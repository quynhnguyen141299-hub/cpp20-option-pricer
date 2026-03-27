#pragma once
/// @file barrier_mc.hpp
/// Monte Carlo engine for barrier options with Brownian bridge correction.
///
/// The Brownian bridge (Beaglehole-Dyer-Jacka) adjusts for the probability
/// that the barrier was crossed between discrete monitoring points:
///   P(hit | S_t, S_{t+dt}) = exp(-2·ln(S_t/B)·ln(S_{t+dt}/B) / (σ²·dt))
///
/// For knock-in options, we use in-out parity:
///   knock_in_price = vanilla_price - knock_out_price
/// This avoids the high-variance "conditional on hitting" estimator.

#include "../core/concepts.hpp"
#include "../core/sobol.hpp"
#include "../core/thread_pool.hpp"
#include "../models/barrier.hpp"
#include "black_scholes.hpp"
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include <format>
#include <random>

namespace pricer {

// ---------------------------------------------------------------------------
// BarrierMCEngine — MC barrier option pricer.
// ---------------------------------------------------------------------------
template <RateModel RMod, VolSurface VMod>
class BarrierMCEngine {
public:
    BarrierMCEngine(RMod rates, VMod vols, MCConfig cfg = {})
        : rates_(std::move(rates)), vols_(std::move(vols)), cfg_(cfg) {}

    BarrierMCEngine(BarrierMCEngine&&) noexcept = default;
    BarrierMCEngine& operator=(BarrierMCEngine&&) noexcept = default;

    /// Price a barrier option.
    [[nodiscard]] Result<PricingResult> price_barrier(
        const BarrierContract& bc, const MarketSnap& m
    ) {
        if (m.S.v <= 0) return std::unexpected(PricingError::at(Err::BadSpot, "S<=0"));
        if (bc.K.v <= 0) return std::unexpected(PricingError::at(Err::BadStrike, "K<=0"));
        if (bc.T.v <= 0) return std::unexpected(PricingError::at(Err::BadExpiry, "T<=0"));
        if (bc.barrier <= 0) return std::unexpected(PricingError::at(Err::BadConfig, "barrier<=0"));

        // For knock-in: use in-out parity
        if (!is_knock_out(bc.barrier_type)) {
            BarrierContract ko_bc = bc;
            ko_bc.barrier_type = is_up_barrier(bc.barrier_type)
                ? BarrierType::UpAndOut : BarrierType::DownAndOut;
            ko_bc.rebate = 0.0;

            auto ko_result = price_barrier(ko_bc, m);
            if (!ko_result) return ko_result;

            auto vanilla = vanilla_price(bc.to_vanilla(), m);
            if (!vanilla) return vanilla;

            return PricingResult{
                .price      = vanilla->price - ko_result->price,
                .std_err    = ko_result->std_err,
                .greeks     = {.delta = vanilla->greeks.delta - ko_result->greeks.delta},
                .elapsed_us = ko_result->elapsed_us + vanilla->elapsed_us,
                .method     = std::format("BarrierMC-InOut {}", str(bc.barrier_type))
            };
        }

        auto t0 = std::chrono::high_resolution_clock::now();

        const double S0    = m.S.v;
        const double K     = bc.K.v;
        const double T     = bc.T.v;
        const double sigma = vols_.iv(K, T);
        const double drift = m.drift();
        const double df_d  = rates_.df(T);
        const double B     = bc.barrier;
        const double rebate = bc.rebate;
        const bool is_up   = is_up_barrier(bc.barrier_type);
        const std::uint32_t steps = std::max(cfg_.steps, 10u);
        const double dt    = T / steps;

        unsigned n_threads = cfg_.n_threads ? cfg_.n_threads
                           : std::max(1u, std::thread::hardware_concurrency());
        std::uint64_t paths_per_thread = cfg_.n_paths / n_threads;

        auto simulate_block = [&](std::uint64_t n_paths, std::uint64_t seed) {
            return simulate_ko_paths(
                S0, K, T, dt, sigma, drift, df_d, B, rebate,
                is_up, bc.type, steps, n_paths, seed
            );
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
            auto total_up = simulate_ko_paths(
                S0 + bump, K, T, dt, sigma, drift, df_d, B, rebate,
                is_up, bc.type, steps,
                std::min(cfg_.n_paths, (std::uint64_t)50000),
                cfg_.seed + 777777);
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
            .method     = std::format("BarrierMC-BB {} {}K paths {}steps",
                                      str(bc.barrier_type), total.count / 1000, steps)
        };
    }

    /// Satisfy Engine concept: price vanilla options (delegates to BS).
    [[nodiscard]] Result<PricingResult> price(
        const Contract& c, const MarketSnap& m
    ) {
        return vanilla_price(c, m);
    }

private:
    [[nodiscard]] Result<PricingResult> vanilla_price(
        const Contract& c, const MarketSnap& m
    ) const {
        if (m.S.v <= 0) return std::unexpected(PricingError::at(Err::BadSpot, "S<=0"));
        if (c.K.v <= 0) return std::unexpected(PricingError::at(Err::BadStrike, "K<=0"));
        if (c.T.v <= 0) return std::unexpected(PricingError::at(Err::BadExpiry, "T<=0"));
        double sigma = vols_.iv(c.K.v, c.T.v);
        if (sigma <= 0) return std::unexpected(PricingError::at(Err::BadVol, "σ<=0"));
        auto gk = garman_kohlhagen(m.S.v, c.K.v, c.T.v, sigma, m.r_d.v, m.r_f.v, c.type);
        return PricingResult{
            .price = gk.price, .std_err = 0.0, .greeks = gk.greeks,
            .elapsed_us = 0.0, .method = "GK-Analytic"
        };
    }

    [[nodiscard]] BatchResult simulate_ko_paths(
        double S0, double K, double T, double dt, double sigma,
        double drift, double df_d, double B, double rebate,
        bool is_up, OptType type, std::uint32_t steps,
        std::uint64_t n_paths, std::uint64_t seed
    ) const noexcept {
        const double sqdt      = std::sqrt(dt);
        const double drift_dt  = (drift - 0.5 * sigma * sigma) * dt;
        const double sigma_sqdt = sigma * sqdt;
        const double sigma2_dt = sigma * sigma * dt;

        std::mt19937_64 rng(seed);
        std::normal_distribution<double> norm(0.0, 1.0);
        std::uniform_real_distribution<double> unif(0.0, 1.0);

        BatchResult br;

        for (std::uint64_t path = 0; path < n_paths; ++path) {
            double S_prev = S0;
            bool knocked = false;

            for (std::uint32_t step = 0; step < steps; ++step) {
                double z = norm(rng);
                double S_next = S_prev * std::exp(drift_dt + sigma_sqdt * z);

                if (!knocked) {
                    // Check barrier at the discrete monitoring point
                    if (is_up && S_next >= B) knocked = true;
                    if (!is_up && S_next <= B) knocked = true;

                    // Brownian bridge continuity correction:
                    // If both S_prev and S_next are on the same side of B,
                    // there's still a probability the path crossed B between them.
                    if (!knocked && sigma2_dt > 1e-15) {
                        double log_prev = std::log(S_prev / B);
                        double log_next = std::log(S_next / B);
                        if (log_prev * log_next > 0) {
                            double p_hit = std::exp(-2.0 * log_prev * log_next / sigma2_dt);
                            if (unif(rng) < p_hit) knocked = true;
                        }
                    }
                }

                S_prev = S_next;
            }

            double payoff;
            if (knocked) {
                payoff = rebate;
            } else {
                payoff = (type == OptType::Call)
                    ? std::max(S_prev - K, 0.0)
                    : std::max(K - S_prev, 0.0);
            }
            double pv = df_d * payoff;

            br.sum    += pv;
            br.sum_sq += pv * pv;
            ++br.count;
        }

        return br;
    }

    RMod rates_;
    VMod vols_;
    MCConfig cfg_;
};

static_assert(Engine<BarrierMCEngine<FlatRate, FlatVol>>);

} // namespace pricer
