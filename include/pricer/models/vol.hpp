#pragma once
/// @file vol.hpp
/// Vol surface models satisfying the VolSurface concept.

#include "../core/concepts.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace pricer {

// ---------------------------------------------------------------------------
// FlatVol — constant σ.
// ---------------------------------------------------------------------------
class FlatVol {
    double s_;
public:
    explicit constexpr FlatVol(Vol s) noexcept : s_(s.v) {}

    [[nodiscard]] constexpr double local_vol(double, double, double) const noexcept { return s_; }
    [[nodiscard]] constexpr double iv(double, double)                const noexcept { return s_; }
};

static_assert(VolSurface<FlatVol>);

// ---------------------------------------------------------------------------
// TermStructureVol — time-dependent ATM vol with variance interpolation.
// Uses total-variance interpolation (market convention) to avoid
// calendar arbitrage: σ²(t)·t must be non-decreasing.
// ---------------------------------------------------------------------------
class TermStructureVol {
public:
    struct Node { double t; double vol; };

    explicit TermStructureVol(std::vector<Node> nodes)
        : nodes_(std::move(nodes))
    {
        std::ranges::sort(nodes_, {}, &Node::t);
    }

    [[nodiscard]] double local_vol(double, double, double t) const { return interp(t); }
    [[nodiscard]] double iv(double, double t)                const { return interp(t); }

private:
    [[nodiscard]] double interp(double t) const {
        if (nodes_.empty()) return 0;
        if (t <= nodes_.front().t) return nodes_.front().vol;
        if (t >= nodes_.back().t)  return nodes_.back().vol;
        auto it = std::ranges::lower_bound(nodes_, t, {}, &Node::t);
        auto prev = std::prev(it);
        // Variance interpolation: total_var = σ² · t
        double v2_prev = prev->vol * prev->vol * prev->t;
        double v2_next = it->vol   * it->vol   * it->t;
        double alpha = (t - prev->t) / (it->t - prev->t);
        double v2_t  = v2_prev + alpha * (v2_next - v2_prev);
        return std::sqrt(v2_t / t);
    }

    std::vector<Node> nodes_;
};

static_assert(VolSurface<TermStructureVol>);

// ---------------------------------------------------------------------------
// StickyStrikeVol — quadratic smile: σ(K,T) = σ_ATM · (1 + β·m + γ·m²)
//   where m = log(K/S_ref).
// Captures risk-reversal (β) and butterfly (γ) in one shot.
// ---------------------------------------------------------------------------
class StickyStrikeVol {
    double atm_, beta_, gamma_, s_ref_;
public:
    StickyStrikeVol(double atm, double rr, double bf, double s_ref)
        : atm_(atm), beta_(rr), gamma_(bf), s_ref_(s_ref) {}

    [[nodiscard]] double local_vol(double S, double K, double) const {
        double m = std::log(K / S);
        return atm_ * (1.0 + beta_ * m + gamma_ * m * m);
    }

    [[nodiscard]] double iv(double K, double) const {
        double m = std::log(K / s_ref_);
        return atm_ * (1.0 + beta_ * m + gamma_ * m * m);
    }
};

static_assert(VolSurface<StickyStrikeVol>);

} // namespace pricer
