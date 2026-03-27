#pragma once
/// @file barrier.hpp
/// Barrier option types: knock-in, knock-out, up, down.

#include "../core/types.hpp"
#include <cstdint>

namespace pricer {

// ---------------------------------------------------------------------------
enum class BarrierType : std::uint8_t {
    UpAndOut,
    DownAndOut,
    UpAndIn,
    DownAndIn
};

[[nodiscard]] constexpr std::string_view str(BarrierType bt) noexcept {
    switch (bt) {
        case BarrierType::UpAndOut:   return "UpAndOut";
        case BarrierType::DownAndOut: return "DownAndOut";
        case BarrierType::UpAndIn:    return "UpAndIn";
        case BarrierType::DownAndIn:  return "DownAndIn";
    }
    return "Unknown";
}

[[nodiscard]] constexpr bool is_knock_out(BarrierType bt) noexcept {
    return bt == BarrierType::UpAndOut || bt == BarrierType::DownAndOut;
}

[[nodiscard]] constexpr bool is_up_barrier(BarrierType bt) noexcept {
    return bt == BarrierType::UpAndOut || bt == BarrierType::UpAndIn;
}

// ---------------------------------------------------------------------------
// BarrierContract — extends Contract with barrier-specific fields.
// ---------------------------------------------------------------------------
struct BarrierContract {
    Strike      K;
    YearFrac    T;
    OptType     type;
    Exercise    exercise     = Exercise::European;
    std::string underlying   = "EURUSD";
    double      barrier      = 0.0;   ///< barrier level
    BarrierType barrier_type = BarrierType::DownAndOut;
    double      rebate       = 0.0;   ///< paid on knock-out

    /// Convert to a vanilla Contract (for in-out parity).
    [[nodiscard]] Contract to_vanilla() const {
        return Contract{K, T, type, exercise, underlying};
    }
};

} // namespace pricer
