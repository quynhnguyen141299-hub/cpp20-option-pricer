#pragma once
/// @file heston_mc.hpp
/// Monte Carlo engine for the Heston stochastic vol model.
///
/// Uses the QE (Quadratic Exponential) discretisation scheme by Andersen (2008)
/// for the variance process.  This is the industry standard because it correctly
/// handles the non-negativity of variance without Euler discretisation artifacts.
///
/// QE scheme summary:
///   Given v_t, compute conditional mean m and variance s² of v_{t+dt}.
///   ψ = s²/m² measures how "spread" the distribution is.
///   - ψ ≤ ψ_crit (≈1.5): use moment-matched quadratic.
///   - ψ > ψ_crit: use exponential approximation.
///
/// Correlation: Cholesky decomposition Z_S = ρ·Z_1 + √(1-ρ²)·Z_2.
/// Delta: bump-and-reprice with dS = 0.01·S.

#include "../core/concepts.hpp"
#include "../core/sobol.hpp"
#include "../core/thread_pool.hpp"
#include "../models/heston.hpp"
#include "black_scholes.hpp"
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include <format>
#include <random>

namespace pricer {

// ---------------------------------------------------------------------------
// HestonMCEngine — Heston MC with QE variance discretisation.
// ---------------------------------------------------------------------------
template <RateModel RMod>
class HestonMCEngine {
public:
    HestonMCEngine(RMod rates, HestonParams params, MCConfig cfg = {})
        : rates_(std::move(rates)), params_(params), cfg_(cfg) {}

    HestonMCEngine(HestonMCEngine&&) noexcept = default;
    HestonMCEngine& operator=(HestonMCEngine&&) noexcept = default;

    [[nodiscard]] Result<PricingResult> price(
        const Contract& c, const MarketSnap& m
    ) {
        if (m.S.v <= 0) return std::unexpected(PricingError::at(Err::BadSpot, "S<=0"));
        if (c.K.v <= 0) return std::unexpected(PricingError::at(Err::BadStrike, "K<=0"));
        if (c.T.v <= 0) return std::unexpected(PricingError::at(Err::BadExpiry, "T<=0"));

        auto t0 = std::chrono::high_resolution_clock::now();

        const double S0    = m.S.v;
        const double K     = c.K.v;
        const double T     = c.T.v;
        const double drift = m.drift();
        const double df_d  = rates_.df(T);
        const std::uint32_t steps = std::max(cfg_.steps, 1u);
        const double dt    = T / steps;

        unsigned n_threads = cfg_.n_threads ? cfg_.n_threads
                           : std::max(1u, std::thread::hardware_concurrency());
        std::uint64_t paths_per_thread = cfg_.n_paths / n_threads;

        auto simulate_block = [&](std::uint64_t n_paths, std::uint64_t seed) {
            return simulate_paths(S0, K, T, dt, drift, df_d, steps,
                                  c.type, n_paths, seed);
        };

        BatchResult total;
        if (n_threads > 1) {
            ThreadPool pool(n_threads);
            std::vector<std::future<BatchResult>> futures;
            futures.reserve(n_threads);
            for (unsigned t = 0; t < n_threads; ++t) {
                std::uint64_t seed = cfg_.seed + t * 100000;
                std::uint64_t n = (t == n_threads - 1)
                    ? (cfg_.n_paths - paths_per_thread * (n_threads - 1))
                    : paths_per_thread;
                futures.push_back(pool.submit(simulate_block, n, seed));
            }
            for (auto& f : futures) total += f.get();
        } else {
            total = simulate_block(cfg_.n_paths, cfg_.seed);
        }

        // Bump-and-reprice delta
        double delta = 0;
        if (total.count > 0) {
            double bump = 0.01 * S0;
            auto total_up = simulate_paths(S0 + bump, K, T, dt, drift, df_d, steps,
                                           c.type, std::min(cfg_.n_paths, (std::uint64_t)50000),
                                           cfg_.seed + 999999);
            delta = (total_up.mean() - total.mean()) / bump;
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        double n = static_cast<double>(total.count);

        Greeks g;
        g.delta = delta;

        return PricingResult{
            .price      = total.mean(),
            .std_err    = std::sqrt(total.variance() / n),
            .greeks     = g,
            .elapsed_us = us,
            .method     = std::format("HestonMC-QE {}K paths {}steps",
                                      total.count / 1000, steps)
        };
    }

private:
    [[nodiscard]] BatchResult simulate_paths(
        double S0, double K, double T, double dt, double drift,
        double df_d, std::uint32_t steps, OptType type,
        std::uint64_t n_paths, std::uint64_t seed
    ) const noexcept {
        const double kappa = params_.kappa;
        const double theta = params_.theta;
        const double xi    = params_.xi;   // vol of vol (σ_v)
        const double rho   = params_.rho;
        const double rho2  = rho * rho;

        // QE precomputed constants for variance
        const double exp_kdt   = std::exp(-kappa * dt);
        const double m_coeff_v = theta * (1.0 - exp_kdt);
        const double s2_c1     = theta * xi * xi * (1.0 - exp_kdt) * (1.0 - exp_kdt) / (2.0 * kappa);
        const double s2_c2     = xi * xi * exp_kdt * (1.0 - exp_kdt) / kappa;
        constexpr double psi_crit = 1.5;

        // Andersen (2008) Section 4: log-spot coefficients for γ₁=γ₂=0.5
        const double gamma1 = 0.5;
        const double gamma2 = 0.5;
        const double K0 = -rho * kappa * theta / xi * dt + drift * dt;
        const double K1 = gamma1 * dt * (kappa * rho / xi - 0.5) - rho / xi;
        const double K2 = gamma2 * dt * (kappa * rho / xi - 0.5) + rho / xi;
        const double K3 = gamma1 * dt * (1.0 - rho2);
        const double K4 = gamma2 * dt * (1.0 - rho2);

        std::mt19937_64 rng(seed);
        std::normal_distribution<double> norm(0.0, 1.0);
        std::uniform_real_distribution<double> unif(0.0, 1.0);

        BatchResult br;

        for (std::uint64_t path = 0; path < n_paths; ++path) {
            double lnS = std::log(S0);
            double v = params_.v0;

            for (std::uint32_t step = 0; step < steps; ++step) {
                double z_v = norm(rng);  // for variance
                double z_s = norm(rng);  // independent, for spot (correlation encoded in K1,K2)

                // QE discretisation for variance
                double m  = v * exp_kdt + m_coeff_v;
                double s2 = v * s2_c2 + s2_c1;
                s2 = std::max(s2, 0.0);
                double psi = (m > 1e-15) ? s2 / (m * m) : 100.0;

                double v_next;
                if (psi <= psi_crit) {
                    double inv_psi = 1.0 / psi;
                    double b2 = 2.0 * inv_psi - 1.0 + std::sqrt(2.0 * inv_psi) * std::sqrt(std::max(2.0 * inv_psi - 1.0, 0.0));
                    double b_qe = std::sqrt(std::max(b2, 0.0));
                    double a_qe = m / (1.0 + b2);
                    v_next = a_qe * (b_qe + z_v) * (b_qe + z_v);
                } else {
                    double p_val = (psi - 1.0) / (psi + 1.0);
                    double beta  = (1.0 - p_val) / std::max(m, 1e-15);
                    double u_exp = unif(rng);
                    if (u_exp <= p_val) {
                        v_next = 0.0;
                    } else {
                        v_next = std::log((1.0 - p_val) / (1.0 - u_exp)) / beta;
                    }
                }
                v_next = std::max(v_next, 0.0);

                // Andersen (2008) log-spot update:
                // ln(S_{t+dt}) = ln(S_t) + K0 + K1·v_t + K2·v_{t+dt} + √(K3·v_t + K4·v_{t+dt})·Z_S
                double var_S = K3 * v + K4 * v_next;
                lnS += K0 + K1 * v + K2 * v_next + std::sqrt(std::max(var_S, 0.0)) * z_s;

                v = v_next;
            }

            double S = std::exp(lnS);

            double payoff = (type == OptType::Call)
                ? std::max(S - K, 0.0)
                : std::max(K - S, 0.0);
            double pv = df_d * payoff;

            br.sum    += pv;
            br.sum_sq += pv * pv;
            ++br.count;
        }

        return br;
    }

    RMod rates_;
    HestonParams params_;
    MCConfig cfg_;
};

static_assert(Engine<HestonMCEngine<FlatRate>>);

} // namespace pricer
