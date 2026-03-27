#pragma once
/// @file types.hpp
/// Domain value types with strong typing to prevent parameter transposition.
/// Every type is trivially copyable and constexpr-constructible.

#include <cstdint>
#include <string>
#include <string_view>
#include <format>
#include <utility>

namespace pricer {

// ---------------------------------------------------------------------------
// StrongDouble — phantom-tagged double that prevents accidental swaps.
// Zero runtime cost: same layout as double, no vtable, no indirection.
// ---------------------------------------------------------------------------
template <typename Tag>
struct StrongDouble {
    double v;
    constexpr explicit StrongDouble(double x) noexcept : v(x) {}
    constexpr operator double() const noexcept { return v; }
    constexpr auto operator<=>(const StrongDouble&) const noexcept = default;
};

using Spot         = StrongDouble<struct SpotTag>;
using Strike       = StrongDouble<struct StrikeTag>;
using Vol          = StrongDouble<struct VolTag>;
using Rate         = StrongDouble<struct RateTag>;
using YearFrac     = StrongDouble<struct YearFracTag>;

// ---------------------------------------------------------------------------
enum class OptType  : std::uint8_t { Call, Put };
enum class Exercise : std::uint8_t { European, Bermudan, American };

[[nodiscard]] constexpr std::string_view str(OptType t) noexcept {
    return t == OptType::Call ? "Call" : "Put";
}

// ---------------------------------------------------------------------------
// Contract — immutable option specification.
// ---------------------------------------------------------------------------
struct Contract {
    Strike    K;
    YearFrac  T;
    OptType   type;
    Exercise  exercise = Exercise::European;
    std::string underlying = "EURUSD";
};

// ---------------------------------------------------------------------------
// MarketSnap — point-in-time market data.  Designed for FX (Garman-Kohlhagen):
//   domestic_rate = ccy2 (USD for EURUSD), foreign_rate = ccy1 (EUR).
//   For equity: foreign_rate = dividend yield.
// ---------------------------------------------------------------------------
struct MarketSnap {
    Spot   S;
    Vol    sigma;
    Rate   r_d;       // domestic / risk-free
    Rate   r_f;       // foreign  / dividend yield

    [[nodiscard]] constexpr double drift() const noexcept {
        return r_d.v - r_f.v;
    }
};

// ---------------------------------------------------------------------------
// Greeks — first- and second-order sensitivities.
// ---------------------------------------------------------------------------
struct Greeks {
    double delta = 0, gamma = 0, vega = 0, theta = 0, rho = 0;

    [[nodiscard]] std::string str() const {
        return std::format("Δ={:+.6f} Γ={:.6f} V={:.4f} Θ={:.4f} ρ={:.4f}",
                           delta, gamma, vega, theta, rho);
    }
};

// ---------------------------------------------------------------------------
// PricingResult — output from any engine.
// ---------------------------------------------------------------------------
struct PricingResult {
    double      price      = 0;
    double      std_err    = 0;   // 0 for analytic
    Greeks      greeks     = {};
    double      elapsed_us = 0;   // microseconds
    std::string method;

    [[nodiscard]] std::string summary() const {
        return std::format("[{}] px={:.6f} se={:.2e} {:.0f}µs\n  {}",
                           method, price, std_err, elapsed_us, greeks.str());
    }
};

// ---------------------------------------------------------------------------
// MCConfig — simulation parameters.
// ---------------------------------------------------------------------------
struct MCConfig {
    std::uint64_t  n_paths         = 100'000;
    std::uint32_t  steps           = 1;       // 1 = terminal only
    std::uint64_t  seed            = 42;
    bool           antithetic      = true;
    bool           control_variate = true;
    double         target_se       = 0.0;     // >0 enables early stop
    std::uint32_t  n_threads       = 0;       // 0 = hardware_concurrency
    std::uint32_t  batch_size      = 10'000;  // paths per coroutine yield
};

} // namespace pricer
