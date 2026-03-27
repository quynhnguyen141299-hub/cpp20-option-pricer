#pragma once
/// @file heston.hpp
/// Heston stochastic volatility model.
///
/// dS = (r_d - r_f) S dt + √v S dW_1
/// dv = κ(θ - v) dt + ξ√v dW_2
/// corr(dW_1, dW_2) = ρ
///
/// Semi-closed-form pricing via Heston (1993) with the Albrecher et al.
/// "little Heston trap" fix, integrated with 64-point Gauss-Legendre.

#include "../core/concepts.hpp"
#include "../engines/black_scholes.hpp"
#include <cmath>
#include <complex>
#include <numbers>

namespace pricer {

// ---------------------------------------------------------------------------
// HestonParams — model parameters.
// ---------------------------------------------------------------------------
struct HestonParams {
    double v0    = 0.04;   ///< initial variance
    double kappa = 1.5;    ///< mean reversion speed
    double theta = 0.04;   ///< long-run variance
    double xi    = 0.3;    ///< vol of vol
    double rho   = -0.7;   ///< correlation W_1, W_2

    /// Feller condition: 2κθ > ξ².  If violated variance can hit zero.
    [[nodiscard]] constexpr bool feller_satisfied() const noexcept {
        return 2.0 * kappa * theta > xi * xi;
    }
};

namespace detail_heston {

// 32-point Gauss-Legendre nodes/weights on [-1,1]
static constexpr int NGL = 32;
static constexpr double gl_x[NGL] = {
    -0.9972638618, -0.9856115115, -0.9647622556, -0.9349060759,
    -0.8963211558, -0.8493676137, -0.7944837960, -0.7321821187,
    -0.6630442669, -0.5877157572, -0.5068999089, -0.4213512761,
    -0.3318686023, -0.2392873623, -0.1444719616, -0.0483076657,
     0.0483076657,  0.1444719616,  0.2392873623,  0.3318686023,
     0.4213512761,  0.5068999089,  0.5877157572,  0.6630442669,
     0.7321821187,  0.7944837960,  0.8493676137,  0.8963211558,
     0.9349060759,  0.9647622556,  0.9856115115,  0.9972638618
};
static constexpr double gl_w[NGL] = {
    0.0070186100,  0.0162743947,  0.0253920653,  0.0342738629,
    0.0428358980,  0.0509980593,  0.0586840935,  0.0658222228,
    0.0723457941,  0.0781938958,  0.0833119242,  0.0876520930,
    0.0911738787,  0.0938443991,  0.0956387201,  0.0965400885,
    0.0965400885,  0.0956387201,  0.0938443991,  0.0911738787,
    0.0876520930,  0.0833119242,  0.0781938958,  0.0723457941,
    0.0658222228,  0.0586840935,  0.0509980593,  0.0428358980,
    0.0342738629,  0.0253920653,  0.0162743947,  0.0070186100
};

/// Heston (1993) integrand for the call price.
/// Uses the formulation from Gatheral (2006), "The Volatility Surface",
/// with the "little Heston trap" (Albrecher et al. 2007) fix.
///
/// Call = S·exp(-r_f·T)·P1 - K·exp(-r_d·T)·P2
/// P_j = 0.5 + (1/π) ∫_0^∞ Re[ e^{-i·φ·ln(K)} · f_j(φ) / (i·φ) ] dφ
///
/// where f_j is the characteristic function for j=1,2 with different parameters.
[[nodiscard]] inline double heston_P(
    int j,  // 1 or 2
    double S, double K, double T, double r_d, double r_f,
    const HestonParams& p
) noexcept {
    using C = std::complex<double>;
    const C i1{0.0, 1.0};

    const double kappa = p.kappa;
    const double theta_val = p.theta;
    const double sigma = p.xi;
    const double rho   = p.rho;
    const double v0    = p.v0;
    const double lnS   = std::log(S);
    const double lnK   = std::log(K);
    const double sigma2 = sigma * sigma;

    // u_j and b_j differ between P1 and P2
    const double u_j = (j == 1) ? 0.5 : -0.5;
    const double b_j = (j == 1) ? (kappa - rho * sigma) : kappa;

    auto integrand = [&](double phi) -> double {
        const C iphi = i1 * phi;
        const C rsi  = C{rho * sigma, 0.0} * iphi;

        // d = sqrt((rho·σ·i·φ - b_j)² - σ²·(2·u_j·i·φ - φ²))
        const C tmp = rsi - C{b_j, 0.0};
        const C d = std::sqrt(tmp * tmp - C{sigma2, 0.0} * (2.0 * C{u_j, 0.0} * iphi - C{phi * phi, 0.0}));

        // "Little Heston trap" fix: use g2 formulation for numerical stability.
        // g = (b_j - rho·σ·i·φ + d) / (b_j - rho·σ·i·φ - d)   [ORIGINAL]
        // Instead use: g = (b_j - rho·σ·i·φ - d) / (b_j - rho·σ·i·φ + d)  [TRAP FIX]
        const C g_num = C{b_j, 0.0} - rsi - d;
        const C g_den = C{b_j, 0.0} - rsi + d;
        const C g = g_num / g_den;

        const C exp_neg_dT = std::exp(-d * C{T, 0.0});

        // D coefficient: with the trap fix
        const C D_val = (C{b_j, 0.0} - rsi - d) / sigma2
            * (1.0 - exp_neg_dT) / (1.0 - g * exp_neg_dT);

        // C coefficient
        const C C_val = C{(r_d - r_f), 0.0} * iphi * C{T, 0.0}
            + C{kappa * theta_val / sigma2, 0.0}
            * ((C{b_j, 0.0} - rsi - d) * C{T, 0.0}
               - 2.0 * std::log((1.0 - g * exp_neg_dT) / (1.0 - g)));

        const C f = std::exp(C_val + D_val * C{v0, 0.0} + iphi * C{lnS, 0.0});

        const C result = std::exp(-iphi * C{lnK, 0.0}) * f / (iphi);
        return result.real();
    };

    // Integrate over [0, u_max] using Gauss-Legendre on multiple subintervals
    auto gl_seg = [&](double a, double b) -> double {
        double half = (b - a) / 2.0;
        double mid  = (b + a) / 2.0;
        double sum  = 0.0;
        for (int k = 0; k < NGL; ++k) {
            double u = half * gl_x[k] + mid;
            if (u < 1e-10) continue;
            sum += gl_w[k] * integrand(u);
        }
        return half * sum;
    };

    double integral = gl_seg(1e-8, 20.0) + gl_seg(20.0, 80.0) + gl_seg(80.0, 200.0);

    double P = 0.5 + integral / std::numbers::pi;
    return std::clamp(P, 0.0, 1.0);
}

/// Heston semi-closed-form call price.
[[nodiscard]] inline double heston_call_price(
    double S, double K, double T, double r_d, double r_f,
    const HestonParams& p
) noexcept {
    double P1 = heston_P(1, S, K, T, r_d, r_f, p);
    double P2 = heston_P(2, S, K, T, r_d, r_f, p);
    double df_d = std::exp(-r_d * T);
    double df_f = std::exp(-r_f * T);
    double call = S * df_f * P1 - K * df_d * P2;
    return std::max(call, 0.0);
}

} // namespace detail_heston

// ---------------------------------------------------------------------------
// HestonVol — satisfies VolSurface concept.
//   local_vol: returns sqrt(v0) as approximation.
//   iv: inverts Heston call price via Newton-Raphson.
// ---------------------------------------------------------------------------
class HestonVol {
    HestonParams params_;
    double r_d_;
    double r_f_;
    double S_;

public:
    HestonVol(HestonParams p, double S, double r_d, double r_f) noexcept
        : params_(p), r_d_(r_d), r_f_(r_f), S_(S) {}

    [[nodiscard]] double local_vol(double, double, double) const noexcept {
        return std::sqrt(params_.v0);
    }

    /// Extract implied vol by pricing with the Heston CF then inverting via Newton.
    [[nodiscard]] double iv(double K, double T) const noexcept {
        double call_px = detail_heston::heston_call_price(S_, K, T, r_d_, r_f_, params_);

        // Newton-Raphson IV extraction
        double sig = std::sqrt(params_.v0);  // initial guess
        for (int i = 0; i < 100; ++i) {
            auto gk = garman_kohlhagen(S_, K, T, sig, r_d_, r_f_, OptType::Call);
            double diff = gk.price - call_px;
            double vega_raw = gk.greeks.vega / 0.01;
            if (std::abs(vega_raw) < 1e-15) break;
            sig -= diff / vega_raw;
            sig = std::clamp(sig, 1e-6, 5.0);
            if (std::abs(diff) < 1e-12) break;
        }
        return sig;
    }

    [[nodiscard]] const HestonParams& params() const noexcept { return params_; }

    /// Direct access to Heston semi-closed-form call price.
    [[nodiscard]] double call_price(double K, double T) const noexcept {
        return detail_heston::heston_call_price(S_, K, T, r_d_, r_f_, params_);
    }

    /// Put price via put-call parity.
    [[nodiscard]] double put_price(double K, double T) const noexcept {
        double c = call_price(K, T);
        double fwd = S_ * std::exp((r_d_ - r_f_) * T);
        return c - (fwd - K) * std::exp(-r_d_ * T);
    }
};

static_assert(VolSurface<HestonVol>);

} // namespace pricer
