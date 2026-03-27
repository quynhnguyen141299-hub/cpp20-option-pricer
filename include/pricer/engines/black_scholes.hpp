#pragma once
/// @file black_scholes.hpp
/// Analytical Garman-Kohlhagen pricer with closed-form Greeks.
/// Also serves as the control variate baseline for the MC engine.

#include "../core/concepts.hpp"
#include <cmath>
#include <numbers>
#include <chrono>

namespace pricer {

namespace detail {
    [[nodiscard]] inline double N(double x)  noexcept { return 0.5 * std::erfc(-x * std::numbers::sqrt2 / 2.0); }
    [[nodiscard]] inline double n(double x)  noexcept { return std::exp(-0.5*x*x) / std::sqrt(2.0 * std::numbers::pi); }
} // namespace detail

/// Closed-form GK price + Greeks.  Used directly and as control variate.
struct GKResult {
    double price;
    Greeks greeks;
};

/// Compute GK price and all first-order Greeks in a single pass.
/// No redundant CDF/PDF evaluations.
[[nodiscard]] inline GKResult garman_kohlhagen(
    double S, double K, double T, double sigma, double r_d, double r_f,
    OptType type
) noexcept {
    const double sqT   = std::sqrt(T);
    const double svT   = sigma * sqT;
    const double d1    = (std::log(S / K) + (r_d - r_f + 0.5 * sigma * sigma) * T) / svT;
    const double d2    = d1 - svT;
    const double df_d  = std::exp(-r_d * T);
    const double df_f  = std::exp(-r_f * T);
    const double Nd1   = detail::N(d1);
    const double Nd2   = detail::N(d2);
    const double nd1   = detail::n(d1);

    GKResult r;
    if (type == OptType::Call) {
        r.price          = S * df_f * Nd1 - K * df_d * Nd2;
        r.greeks.delta   = df_f * Nd1;
        r.greeks.gamma   = df_f * nd1 / (S * svT);
        r.greeks.vega    = S * df_f * nd1 * sqT * 0.01;          // per 1% vol
        r.greeks.theta   = (-(S * df_f * nd1 * sigma) / (2.0 * sqT)
                            - r_d * K * df_d * Nd2
                            + r_f * S * df_f * Nd1) / 365.0;     // per day
        r.greeks.rho     = K * T * df_d * Nd2 * 0.01;            // per 1% rate
    } else {
        double Nmd1 = detail::N(-d1);
        double Nmd2 = detail::N(-d2);
        r.price          = K * df_d * Nmd2 - S * df_f * Nmd1;
        r.greeks.delta   = -df_f * Nmd1;
        r.greeks.gamma   = df_f * nd1 / (S * svT);
        r.greeks.vega    = S * df_f * nd1 * sqT * 0.01;
        r.greeks.theta   = (-(S * df_f * nd1 * sigma) / (2.0 * sqT)
                            + r_d * K * df_d * Nmd2
                            - r_f * S * df_f * Nmd1) / 365.0;
        r.greeks.rho     = -K * T * df_d * Nmd2 * 0.01;
    }
    return r;
}

// ---------------------------------------------------------------------------
// BSEngine — wraps garman_kohlhagen into the Engine concept.
// ---------------------------------------------------------------------------
class BSEngine {
public:
    [[nodiscard]] Result<PricingResult> price(
        const Contract& c, const MarketSnap& m
    ) const {
        auto v = validate(c, m);
        if (!v) return std::unexpected(v.error());

        auto t0 = std::chrono::high_resolution_clock::now();

        auto gk = garman_kohlhagen(m.S.v, c.K.v, c.T.v, m.sigma.v,
                                    m.r_d.v, m.r_f.v, c.type);

        auto t1 = std::chrono::high_resolution_clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();

        return PricingResult{
            .price      = gk.price,
            .std_err    = 0.0,
            .greeks     = gk.greeks,
            .elapsed_us = us,
            .method     = "GK-Analytic"
        };
    }

    /// Newton-Raphson implied vol solver.
    [[nodiscard]] Result<double> implied_vol(
        const Contract& c, const MarketSnap& m, double target_px,
        double tol = 1e-10, int max_iter = 100
    ) const {
        if (target_px <= 0)
            return std::unexpected(PricingError::at(Err::BadConfig, "target_px <= 0"));

        double sig = 0.20;
        for (int i = 0; i < max_iter; ++i) {
            auto gk = garman_kohlhagen(m.S.v, c.K.v, c.T.v, sig,
                                        m.r_d.v, m.r_f.v, c.type);
            double diff = gk.price - target_px;
            double vega_raw = gk.greeks.vega / 0.01;
            if (std::abs(vega_raw) < 1e-15)
                return std::unexpected(PricingError::at(Err::NoConverge, "vega~0"));
            sig -= diff / vega_raw;
            sig = std::clamp(sig, 1e-6, 5.0);
            if (std::abs(diff) < tol) return sig;
        }
        return std::unexpected(PricingError::at(Err::NoConverge, "IV max_iter"));
    }

private:
    [[nodiscard]] static Result<void> validate(const Contract& c, const MarketSnap& m) {
        if (auto r = require(m.S.v > 0,     Err::BadSpot,   "S<=0"); !r) return r;
        if (auto r = require(c.K.v > 0,     Err::BadStrike, "K<=0"); !r) return r;
        if (auto r = require(m.sigma.v > 0, Err::BadVol,    "σ<=0"); !r) return r;
        if (auto r = require(c.T.v > 0,     Err::BadExpiry, "T<=0"); !r) return r;
        return {};
    }
};

static_assert(Engine<BSEngine>);

} // namespace pricer
