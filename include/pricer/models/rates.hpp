#pragma once
/// @file rates.hpp
/// Rate models satisfying the RateModel concept.

#include "../core/concepts.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace pricer {

// ---------------------------------------------------------------------------
// FlatRate — constant continuously-compounded rate.
// Typical use: unit tests, quick sanity checks.
// ---------------------------------------------------------------------------
class FlatRate {
    double r_;
public:
    explicit constexpr FlatRate(Rate r) noexcept : r_(r.v) {}

    [[nodiscard]] constexpr double df(double t)  const noexcept { return std::exp(-r_ * t); }
    [[nodiscard]] constexpr double fwd(double)   const noexcept { return r_; }
};

static_assert(RateModel<FlatRate>);

// ---------------------------------------------------------------------------
// PiecewiseCurve — linearly-interpolated zero curve.
// Constructed from bootstrapped nodes (e.g. OIS/SOFR strip).
// Owns its node vector (RAII), supports move for efficient transfer.
// ---------------------------------------------------------------------------
class PiecewiseCurve {
public:
    struct Node { double t; double r; };

    explicit PiecewiseCurve(std::vector<Node> nodes)
        : nodes_(std::move(nodes))
    {
        std::ranges::sort(nodes_, {}, &Node::t);
    }

    PiecewiseCurve(PiecewiseCurve&&) noexcept = default;
    PiecewiseCurve& operator=(PiecewiseCurve&&) noexcept = default;
    PiecewiseCurve(const PiecewiseCurve&) = default;
    PiecewiseCurve& operator=(const PiecewiseCurve&) = default;

    [[nodiscard]] double zero_rate(double t) const noexcept {
        if (nodes_.empty()) return 0;
        if (t <= nodes_.front().t) return nodes_.front().r;
        if (t >= nodes_.back().t)  return nodes_.back().r;
        auto it = std::ranges::lower_bound(nodes_, t, {}, &Node::t);
        auto prev = std::prev(it);
        double alpha = (t - prev->t) / (it->t - prev->t);
        return prev->r + alpha * (it->r - prev->r);
    }

    [[nodiscard]] double df(double t) const noexcept {
        return std::exp(-zero_rate(t) * t);
    }

    [[nodiscard]] double fwd(double t) const noexcept {
        constexpr double dt = 1.0 / 365.0;
        return -(std::log(df(t + dt)) - std::log(df(t))) / dt;
    }

private:
    std::vector<Node> nodes_;
};

static_assert(RateModel<PiecewiseCurve>);

// ---------------------------------------------------------------------------
// FXRatePair — pairs domestic + foreign rate models for GK drift.
// Templated: works with any combination of RateModel implementations.
// ---------------------------------------------------------------------------
template <RateModel D, RateModel F>
class FXRatePair {
    D dom_;
    F fgn_;
public:
    FXRatePair(D d, F f) : dom_(std::move(d)), fgn_(std::move(f)) {}

    [[nodiscard]] double df(double t)       const { return dom_.df(t); }
    [[nodiscard]] double fwd(double t)      const { return dom_.fwd(t); }
    [[nodiscard]] double fx_drift(double t) const { return dom_.fwd(t) - fgn_.fwd(t); }

    [[nodiscard]] const D& dom() const noexcept { return dom_; }
    [[nodiscard]] const F& fgn() const noexcept { return fgn_; }
};

static_assert(RateModel<FXRatePair<FlatRate, FlatRate>>);

} // namespace pricer
