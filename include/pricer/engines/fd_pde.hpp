#pragma once
/// @file fd_pde.hpp
/// 1D Finite-Difference PDE solver for the Black-Scholes equation.
///
/// Features:
///   - Non-uniform sinh grid concentrated near the strike.
///   - θ-scheme time stepping (θ=0 explicit, θ=0.5 Crank-Nicolson, θ=1 implicit).
///   - Thomas algorithm for tridiagonal solve — the 1D specialisation of
///     Kronecker product structure that generalises to 2D ADI schemes.
///   - American option support via early exercise constraint (max operator).
///   - Greeks extraction directly from the grid.
///
/// The Kronecker connection: in 2D (e.g. Heston PDE for S and v), the finite
/// difference discretisation yields A ⊗ I + I ⊗ B + cross terms.  For 1D BS,
/// this reduces to a single tridiagonal system solvable by the Thomas algorithm.
/// The infrastructure here is the building block for that 2D extension.

#include "../core/concepts.hpp"
#include "black_scholes.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <format>

namespace pricer {

// ---------------------------------------------------------------------------
// FDConfig — finite difference grid and solver parameters.
// ---------------------------------------------------------------------------
struct FDConfig {
    std::uint32_t n_space     = 200;   ///< grid points in S direction
    std::uint32_t n_time      = 200;   ///< time steps
    double        s_min_mult  = 0.0;   ///< S_min = s_min_mult * S (0 recommended)
    double        s_max_mult  = 3.0;   ///< S_max = s_max_mult * S
    double        theta       = 0.5;   ///< 0=explicit, 0.5=Crank-Nicolson, 1=implicit
};

namespace detail_fd {

/// Thomas algorithm for tridiagonal system.
/// a[1..n-1] sub-diagonal, b[0..n-1] diagonal, c[0..n-2] super-diagonal.
/// d[0..n-1] = RHS on input, solution on output.
inline void thomas_solve(
    const std::vector<double>& a,
    std::vector<double>& b,
    const std::vector<double>& c,
    std::vector<double>& d
) noexcept {
    int n = static_cast<int>(b.size());
    // Forward sweep
    for (int i = 1; i < n; ++i) {
        double w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    // Back substitution
    d[n - 1] /= b[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
    }
}

/// Build a non-uniform sinh grid concentrated near the strike K.
inline std::vector<double> build_grid(
    double S_min, double S_max, double K, std::uint32_t N
) {
    std::vector<double> S(N);

    // Use a sinh stretching to concentrate points near K.
    // Map uniform ξ ∈ [0,1] → S via:
    //   S(ξ) = A + B·sinh(C·(ξ - ξ_K))
    // where ξ_K = (K - S_min)/(S_max - S_min) centers the concentration at K.
    double range = S_max - S_min;
    double xi_K = std::clamp((K - S_min) / range, 0.1, 0.9);
    double alpha = 3.0;  // concentration parameter

    // Compute mapping constants
    double sinh_low  = std::sinh(-alpha * xi_K);
    double sinh_high = std::sinh(alpha * (1.0 - xi_K));
    double B = range / (sinh_high - sinh_low);
    double A = S_min - B * sinh_low;

    for (std::uint32_t i = 0; i < N; ++i) {
        double xi = static_cast<double>(i) / static_cast<double>(N - 1);
        S[i] = A + B * std::sinh(alpha * (xi - xi_K));
    }
    S[0] = S_min;
    S[N - 1] = S_max;

    return S;
}

} // namespace detail_fd

// ---------------------------------------------------------------------------
// FDEngine — finite difference PDE solver satisfying Engine concept.
// ---------------------------------------------------------------------------
class FDEngine {
    FDConfig cfg_;

public:
    explicit FDEngine(FDConfig cfg = {}) : cfg_(cfg) {}

    FDEngine(FDEngine&&) noexcept = default;
    FDEngine& operator=(FDEngine&&) noexcept = default;

    [[nodiscard]] Result<PricingResult> price(
        const Contract& c, const MarketSnap& m
    ) const {
        if (m.S.v <= 0) return std::unexpected(PricingError::at(Err::BadSpot, "S<=0"));
        if (c.K.v <= 0) return std::unexpected(PricingError::at(Err::BadStrike, "K<=0"));
        if (m.sigma.v <= 0) return std::unexpected(PricingError::at(Err::BadVol, "σ<=0"));
        if (c.T.v <= 0) return std::unexpected(PricingError::at(Err::BadExpiry, "T<=0"));

        auto t0 = std::chrono::high_resolution_clock::now();

        const double S0    = m.S.v;
        const double K     = c.K.v;
        const double T     = c.T.v;
        const double sigma = m.sigma.v;
        const double r     = m.r_d.v;
        const double q     = m.r_f.v;
        const bool is_call = (c.type == OptType::Call);
        const bool is_american = (c.exercise == Exercise::American);

        const int N = static_cast<int>(cfg_.n_space);
        const int M = static_cast<int>(cfg_.n_time);
        const double dt = T / M;
        const double th = cfg_.theta;

        // Build spatial grid
        double S_min = cfg_.s_min_mult * S0;
        double S_max = cfg_.s_max_mult * S0;
        if (S_min < 0) S_min = 0;

        auto S = detail_fd::build_grid(S_min, S_max, K, N);

        // Terminal condition: payoff at expiry
        std::vector<double> V(N);
        for (int i = 0; i < N; ++i) {
            V[i] = is_call ? std::max(S[i] - K, 0.0) : std::max(K - S[i], 0.0);
        }

        // Pre-compute the PDE coefficients at each interior grid point.
        // BS PDE (backward): ∂V/∂t + 0.5σ²S²·∂²V/∂S² + (r-q)S·∂V/∂S - rV = 0
        // Discretise second derivative and first derivative on the non-uniform grid.
        //
        // At grid point i with S[i-1], S[i], S[i+1]:
        //   h_m = S[i] - S[i-1],  h_p = S[i+1] - S[i]
        //   ∂V/∂S ≈ (h_m²·V[i+1] - h_p²·V[i-1] + (h_p²-h_m²)·V[i]) / (h_m·h_p·(h_m+h_p))
        //   ∂²V/∂S² ≈ 2·(h_m·V[i+1] + h_p·V[i-1] - (h_m+h_p)·V[i]) / (h_m·h_p·(h_m+h_p))
        //
        // LV[i] = alpha_i·V[i-1] + beta_i·V[i] + gamma_i·V[i+1]
        std::vector<double> al(N, 0.0), be(N, 0.0), ga(N, 0.0);
        for (int i = 1; i < N - 1; ++i) {
            double h_m = S[i] - S[i - 1];
            double h_p = S[i + 1] - S[i];
            double h_sum = h_m + h_p;

            double coeff_diff = 0.5 * sigma * sigma * S[i] * S[i];
            double coeff_conv = (r - q) * S[i];

            // Second derivative coefficients
            double d2_lower = 2.0 * coeff_diff / (h_m * h_sum);
            double d2_upper = 2.0 * coeff_diff / (h_p * h_sum);
            double d2_mid   = -(d2_lower + d2_upper);

            // First derivative coefficients (central difference for non-uniform grid)
            double d1_lower = -coeff_conv * h_p / (h_m * h_sum);
            double d1_upper = coeff_conv * h_m / (h_p * h_sum);
            double d1_mid   = coeff_conv * (h_p - h_m) / (h_m * h_p);

            al[i] = d2_lower + d1_lower;
            ga[i] = d2_upper + d1_upper;
            be[i] = d2_mid + d1_mid - r;
        }

        // Time stepping backward from T to 0.
        // At each step:  (I - th·dt·L)·V^new = (I + (1-th)·dt·L)·V^old
        // where L is the spatial operator.
        for (int step = 0; step < M; ++step) {
            double tau = (step + 1) * dt;  // time to maturity remaining

            // Build RHS: (I + (1-th)·dt·L)·V
            std::vector<double> rhs(N);
            for (int i = 1; i < N - 1; ++i) {
                rhs[i] = V[i] + (1.0 - th) * dt * (al[i] * V[i-1] + be[i] * V[i] + ga[i] * V[i+1]);
            }

            // Boundary conditions
            if (is_call) {
                rhs[0]     = 0.0;
                rhs[N - 1] = S[N-1] - K * std::exp(-r * (T - tau));
            } else {
                rhs[0]     = K * std::exp(-r * (T - tau)) - S[0];
                if (rhs[0] < 0) rhs[0] = 0;
                rhs[N - 1] = 0.0;
            }

            // Build LHS tridiagonal: (I - th·dt·L)
            std::vector<double> a_tri(N, 0.0);  // sub-diag
            std::vector<double> b_tri(N, 1.0);  // diag
            std::vector<double> c_tri(N, 0.0);  // super-diag

            for (int i = 1; i < N - 1; ++i) {
                a_tri[i] = -th * dt * al[i];
                b_tri[i] = 1.0 - th * dt * be[i];
                c_tri[i] = -th * dt * ga[i];
            }

            // Solve
            detail_fd::thomas_solve(a_tri, b_tri, c_tri, rhs);
            V = rhs;

            // American: enforce early exercise constraint
            if (is_american) {
                for (int i = 0; i < N; ++i) {
                    double intrinsic = is_call
                        ? std::max(S[i] - K, 0.0)
                        : std::max(K - S[i], 0.0);
                    V[i] = std::max(V[i], intrinsic);
                }
            }
        }

        // Interpolate to find V(S0) and extract Greeks.
        // Find grid indices bracketing S0.
        int idx = 1;
        for (int i = 1; i < N; ++i) {
            if (S[i] >= S0) { idx = i; break; }
        }
        if (idx < 2) idx = 2;
        if (idx >= N - 2) idx = N - 3;

        // Quadratic interpolation using three points around S0
        // for better accuracy.
        int i0 = idx - 1, i1 = idx, i2 = idx + 1;
        double S0_val = S0;

        // Lagrange basis
        double L0 = ((S0_val - S[i1]) * (S0_val - S[i2])) / ((S[i0] - S[i1]) * (S[i0] - S[i2]));
        double L1 = ((S0_val - S[i0]) * (S0_val - S[i2])) / ((S[i1] - S[i0]) * (S[i1] - S[i2]));
        double L2 = ((S0_val - S[i0]) * (S0_val - S[i1])) / ((S[i2] - S[i0]) * (S[i2] - S[i1]));
        double price_fd = V[i0] * L0 + V[i1] * L1 + V[i2] * L2;

        // Extract Greeks from the grid near idx
        Greeks greeks;
        if (idx >= 2 && idx < N - 2) {
            // Delta: central difference
            greeks.delta = (V[idx + 1] - V[idx - 1]) / (S[idx + 1] - S[idx - 1]);

            // Gamma: second derivative
            double h_m = S[idx] - S[idx - 1];
            double h_p = S[idx + 1] - S[idx];
            greeks.gamma = 2.0 * (h_m * V[idx + 1] + h_p * V[idx - 1] - (h_m + h_p) * V[idx])
                         / (h_m * h_p * (h_m + h_p));

            // Theta from PDE: θ = rV - 0.5σ²S²Γ - (r-q)SΔ
            greeks.theta = -(0.5 * sigma * sigma * S[idx] * S[idx] * greeks.gamma
                           + (r - q) * S[idx] * greeks.delta
                           - r * V[idx]) / 365.0;
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();

        return PricingResult{
            .price      = price_fd,
            .std_err    = 0.0,
            .greeks     = greeks,
            .elapsed_us = us,
            .method     = std::format("FD-θ{:.1f} {}×{} {}",
                cfg_.theta, N, M, is_american ? "American" : "European")
        };
    }
};

static_assert(Engine<FDEngine>);

} // namespace pricer
