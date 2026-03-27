#pragma once
/// @file signal.hpp
/// Alpha signal framework for execution timing.
///
/// In practice, execution algos don't just slice mechanically — they
/// incorporate short-term signals to decide whether to be more or less
/// aggressive.  Examples:
///   - Spread signal: trade when spread narrows
///   - Momentum signal: accelerate if price is moving against us
///   - Mean-reversion signal: pause if price will revert
///
/// Design: The Signal concept constrains any predictor that maps a
/// MarketTick to a [-1, +1] score.  Algos can query a signal and
/// adjust slice size accordingly.

#include "order_book.hpp"
#include <concepts>
#include <cmath>
#include <deque>
#include <numeric>

namespace pricer::execution {

// ---------------------------------------------------------------------------
// Signal concept — any short-term alpha predictor.
// ---------------------------------------------------------------------------
template <typename S>
concept Signal = requires(S& s, const MarketTick& tick) {
    { s.score(tick) } -> std::convertible_to<double>;  // in [-1, +1]
    { s.name() }      -> std::convertible_to<std::string>;
};

// ---------------------------------------------------------------------------
// SpreadSignal — trade when spread is tighter than average.
//
// Score > 0 when current spread < rolling average → favorable to trade.
// Score < 0 when spread is wide → wait.
// ---------------------------------------------------------------------------
class SpreadSignal {
public:
    explicit SpreadSignal(int lookback = 100) : lookback_(lookback) {}

    [[nodiscard]] double score(const MarketTick& tick) {
        spreads_.push_back(tick.spread_bps());
        if (static_cast<int>(spreads_.size()) > lookback_)
            spreads_.pop_front();

        if (spreads_.size() < 2) return 0;

        double avg = std::accumulate(spreads_.begin(), spreads_.end(), 0.0)
                   / static_cast<double>(spreads_.size());

        if (avg < 1e-10) return 0;
        // Normalise: (avg - current) / avg, clamped to [-1, 1]
        double raw = (avg - tick.spread_bps()) / avg;
        return std::clamp(raw, -1.0, 1.0);
    }

    [[nodiscard]] std::string name() const { return "SpreadSignal"; }

private:
    int lookback_;
    std::deque<double> spreads_;
};

static_assert(Signal<SpreadSignal>);

// ---------------------------------------------------------------------------
// MomentumSignal — detect short-term price trend.
//
// Score > 0 when price is moving in our favor (or against us if buying).
// Uses simple rolling return over lookback window.
// ---------------------------------------------------------------------------
class MomentumSignal {
public:
    explicit MomentumSignal(int lookback = 50) : lookback_(lookback) {}

    [[nodiscard]] double score(const MarketTick& tick) {
        mids_.push_back(tick.mid);
        if (static_cast<int>(mids_.size()) > lookback_)
            mids_.pop_front();

        if (mids_.size() < 2) return 0;

        double ret = (mids_.back() - mids_.front()) / mids_.front();
        // Normalise to [-1, 1] assuming daily returns rarely exceed 1%
        return std::clamp(ret / 0.01, -1.0, 1.0);
    }

    [[nodiscard]] std::string name() const { return "MomentumSignal"; }

private:
    int lookback_;
    std::deque<double> mids_;
};

static_assert(Signal<MomentumSignal>);

// ---------------------------------------------------------------------------
// CompositeSignal — weighted combination of signals.
//
// Demonstrates generic programming: works with any mix of Signal types
// via type erasure (std::function), or at compile time via templates.
// ---------------------------------------------------------------------------
class CompositeSignal {
public:
    using SignalFn = std::function<double(const MarketTick&)>;

    void add(std::string name, double weight, SignalFn fn) {
        components_.push_back({std::move(name), weight, std::move(fn)});
    }

    [[nodiscard]] double score(const MarketTick& tick) {
        double total = 0, weight_sum = 0;
        for (auto& [n, w, fn] : components_) {
            total += w * fn(tick);
            weight_sum += std::abs(w);
        }
        return weight_sum > 0 ? std::clamp(total / weight_sum, -1.0, 1.0) : 0;
    }

    [[nodiscard]] std::string name() const { return "Composite"; }

private:
    struct Component {
        std::string name;
        double weight;
        SignalFn fn;
    };
    std::vector<Component> components_;
};

static_assert(Signal<CompositeSignal>);

} // namespace pricer::execution
