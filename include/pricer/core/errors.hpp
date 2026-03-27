#pragma once
/// @file errors.hpp
/// Monadic error handling via std::expected.
/// Every fallible function returns Result<T>.  No exceptions on the hot path.

#include <expected>
#include <string>
#include <format>
#include <source_location>

namespace pricer {

enum class Err : int {
    BadSpot      = 1001,
    BadStrike    = 1002,
    BadVol       = 1003,
    BadExpiry    = 1004,
    BadConfig    = 1005,
    Overflow     = 2001,
    NoConverge   = 2002,
};

struct PricingError {
    Err         code;
    std::string msg;
    std::string loc;   // file:line auto-captured

    static PricingError at(Err c, std::string m,
                           std::source_location s = std::source_location::current()) {
        return {c, std::move(m),
                std::format("{}:{}", s.file_name(), s.line())};
    }

    [[nodiscard]] std::string str() const {
        return std::format("[E{}] {} ({})", static_cast<int>(code), msg, loc);
    }
};

template <typename T>
using Result = std::expected<T, PricingError>;

// Validate a precondition; returns unexpected on failure.
inline auto require(bool cond, Err code, std::string msg,
                    std::source_location s = std::source_location::current())
    -> Result<void>
{
    if (!cond) return std::unexpected(PricingError::at(code, std::move(msg), s));
    return {};
}

} // namespace pricer
