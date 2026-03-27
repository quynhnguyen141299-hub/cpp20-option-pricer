#pragma once
/// @file sobol.hpp
/// Sobol quasi-random sequence generator (dimension-1).
///
/// Why Sobol over Mersenne Twister for MC pricing:
///   Sobol is a low-discrepancy sequence: it fills the unit hypercube more
///   uniformly than pseudo-random, giving O(1/N · (log N)^d) error instead
///   of O(1/√N).  For d=1 terminal-value pricing this means ~10x fewer paths
///   for the same accuracy.
///
/// Implementation: Gray-code optimization of the Sobol direction numbers.
/// Uses Joe-Kuo direction numbers for the first dimension.
/// RAII: owns its state, move-only (no accidental sequence correlation).

#include <cstdint>
#include <array>
#include <cmath>
#include <bit>

namespace pricer {

class SobolEngine {
public:
    static constexpr int MAX_BITS = 52;  // double mantissa precision

    explicit SobolEngine(std::uint64_t skip = 0) noexcept {
        init_direction_numbers();
        // Skip initial points (scrambling)
        for (std::uint64_t i = 0; i < skip; ++i) (void)next();
    }

    // Move-only: avoid accidental sequence duplication
    SobolEngine(SobolEngine&&) noexcept = default;
    SobolEngine& operator=(SobolEngine&&) noexcept = default;
    SobolEngine(const SobolEngine&) = delete;
    SobolEngine& operator=(const SobolEngine&) = delete;

    /// Return next point in [0, 1) and advance.
    [[nodiscard]] double next() noexcept {
        // Gray code: find position of rightmost zero bit in index_
        unsigned c = std::countr_one(index_);
        if (c >= MAX_BITS) c = 0;  // overflow guard
        state_ ^= directions_[c];
        ++index_;
        return static_cast<double>(state_) / static_cast<double>(1ULL << MAX_BITS);
    }

    /// Convert uniform [0,1) to standard normal via rational approximation
    /// (Beasley-Springer-Moro inverse normal).
    [[nodiscard]] static double inv_normal(double u) noexcept {
        // Beasley-Springer-Moro algorithm
        static constexpr double a[] = {
            2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637
        };
        static constexpr double b[] = {
            -8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833
        };
        static constexpr double c[] = {
            0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
            0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
            0.0000321767881768, 0.0000002888167364, 0.0000003960315187
        };

        double y = u - 0.5;
        if (std::abs(y) < 0.42) {
            double r = y * y;
            double x = y * (((a[3]*r + a[2])*r + a[1])*r + a[0])
                         / ((((b[3]*r + b[2])*r + b[1])*r + b[0])*r + 1.0);
            return x;
        }

        double r = (y < 0.0) ? u : (1.0 - u);
        double s = std::log(-std::log(r));
        double x = c[0] + s*(c[1] + s*(c[2] + s*(c[3] + s*(c[4]
                 + s*(c[5] + s*(c[6] + s*(c[7] + s*c[8])))))));
        return (y < 0.0) ? -x : x;
    }

    void reset() noexcept { state_ = 0; index_ = 0; }

private:
    void init_direction_numbers() noexcept {
        // Joe-Kuo direction numbers for dimension 1
        // v[i] = 2^(MAX_BITS - i - 1) for the trivial first dimension
        for (int i = 0; i < MAX_BITS; ++i) {
            directions_[i] = 1ULL << (MAX_BITS - 1 - i);
        }
    }

    std::array<std::uint64_t, MAX_BITS> directions_{};
    std::uint64_t state_ = 0;
    std::uint64_t index_ = 0;
};

} // namespace pricer
