#pragma once
/// @file algo.hpp
/// Execution algorithms: TWAP, VWAP, Implementation Shortfall (IS).
///
/// Design: Each algo satisfies the ExecutionAlgo concept — it receives
/// a MarketTick and returns an optional order slice (quantity to execute
/// now).  The backtest harness drives the clock; the algo decides when
/// and how much to trade.
///
/// This mirrors how real algo execution works:
///   - TWAP: split the order evenly over time (minimise timing risk)
///   - VWAP: weight slices by expected volume profile (track the benchmark)
///   - IS:   front-load to minimise arrival-price shortfall (Almgren-Chriss)
///
/// All algos are RAII, move-only, and stateful (track residual quantity).

#include "order_book.hpp"
#include "../core/concepts.hpp"
#include <optional>
#include <cmath>
#include <vector>
#include <numeric>
#include <format>
#include <algorithm>

namespace pricer::execution {

// ---------------------------------------------------------------------------
// ExecutionOrder — what the client wants executed.
// ---------------------------------------------------------------------------
struct ExecutionOrder {
    double  total_qty     = 10e6;     // e.g. 10M EUR
    Side    side          = Side::Ask; // buying = Ask
    double  start_time_us = 0;
    double  end_time_us   = 3600e6;   // 1 hour in microseconds
    std::string label     = "EURUSD";
};

// ---------------------------------------------------------------------------
// SliceDecision — algo's output at each tick.
// ---------------------------------------------------------------------------
struct SliceDecision {
    double quantity  = 0;      // how much to trade this tick (0 = skip)
    bool   is_final = false;   // last slice?
};

// ---------------------------------------------------------------------------
// ExecutionAlgo concept — any algo must provide:
//   - on_tick(MarketTick) → optional<SliceDecision>
//   - filled() → total quantity filled so far
//   - name() → string identifier
// ---------------------------------------------------------------------------
template <typename A>
concept ExecutionAlgo = requires(A& a, const MarketTick& tick) {
    { a.on_tick(tick) }  -> std::same_as<std::optional<SliceDecision>>;
    { a.filled() }       -> std::convertible_to<double>;
    { a.name() }         -> std::convertible_to<std::string>;
};

// ---------------------------------------------------------------------------
// TWAP — Time-Weighted Average Price.
//
// Splits the parent order into equal-sized slices at regular intervals.
// Simple, predictable, widely used for low-urgency orders.
// ---------------------------------------------------------------------------
class TWAPAlgo {
public:
    explicit TWAPAlgo(ExecutionOrder order, int n_slices = 20)
        : order_(std::move(order))
        , n_slices_(n_slices)
        , slice_qty_(order_.total_qty / n_slices)
        , slice_interval_us_((order_.end_time_us - order_.start_time_us) / n_slices)
        , next_slice_time_(order_.start_time_us + slice_interval_us_)
    {}

    TWAPAlgo(TWAPAlgo&&) noexcept = default;
    TWAPAlgo& operator=(TWAPAlgo&&) noexcept = default;

    [[nodiscard]] std::optional<SliceDecision> on_tick(const MarketTick& tick) {
        if (filled_ >= order_.total_qty) return std::nullopt;
        if (tick.timestamp_us < next_slice_time_) return std::nullopt;

        next_slice_time_ += slice_interval_us_;
        double qty = std::min(slice_qty_, order_.total_qty - filled_);
        filled_ += qty;

        return SliceDecision{qty, filled_ >= order_.total_qty};
    }

    [[nodiscard]] double filled() const noexcept { return filled_; }
    [[nodiscard]] std::string name() const { return "TWAP"; }

private:
    ExecutionOrder order_;
    int    n_slices_;
    double slice_qty_;
    double slice_interval_us_;
    double next_slice_time_;
    double filled_ = 0;
};

static_assert(ExecutionAlgo<TWAPAlgo>);

// ---------------------------------------------------------------------------
// VWAP — Volume-Weighted Average Price.
//
// Weights slices by an expected intraday volume profile (U-shaped curve:
// heavy at open and close, light midday).  Tracks the VWAP benchmark
// by participating proportionally to market volume.
// ---------------------------------------------------------------------------
class VWAPAlgo {
public:
    /// volume_profile: fraction of daily volume in each bucket (must sum to ~1).
    /// If empty, defaults to a synthetic U-shaped profile.
    explicit VWAPAlgo(ExecutionOrder order, int n_buckets = 20,
                      std::vector<double> volume_profile = {})
        : order_(std::move(order)), n_buckets_(n_buckets)
    {
        if (volume_profile.empty()) {
            // Synthetic U-shape: heavier at start and end
            volume_profile.resize(n_buckets);
            for (int i = 0; i < n_buckets; ++i) {
                double t = static_cast<double>(i) / (n_buckets - 1);  // 0 to 1
                volume_profile[i] = 1.0 + 2.0 * (t - 0.5) * (t - 0.5);  // U-shape
            }
            double sum = std::accumulate(volume_profile.begin(), volume_profile.end(), 0.0);
            for (auto& v : volume_profile) v /= sum;
        }

        profile_ = std::move(volume_profile);
        bucket_duration_us_ = (order_.end_time_us - order_.start_time_us) / n_buckets;
    }

    VWAPAlgo(VWAPAlgo&&) noexcept = default;
    VWAPAlgo& operator=(VWAPAlgo&&) noexcept = default;

    [[nodiscard]] std::optional<SliceDecision> on_tick(const MarketTick& tick) {
        if (filled_ >= order_.total_qty) return std::nullopt;

        int bucket = static_cast<int>(
            (tick.timestamp_us - order_.start_time_us) / bucket_duration_us_);
        bucket = std::clamp(bucket, 0, n_buckets_ - 1);

        if (bucket <= last_bucket_) return std::nullopt;  // already traded this bucket
        last_bucket_ = bucket;

        double target_qty = profile_[bucket] * order_.total_qty;
        target_qty = std::min(target_qty, order_.total_qty - filled_);
        if (target_qty < 1.0) return std::nullopt;

        filled_ += target_qty;
        return SliceDecision{target_qty, filled_ >= order_.total_qty};
    }

    [[nodiscard]] double filled() const noexcept { return filled_; }
    [[nodiscard]] std::string name() const { return "VWAP"; }

private:
    ExecutionOrder order_;
    int    n_buckets_;
    std::vector<double> profile_;
    double bucket_duration_us_ = 0;
    int    last_bucket_         = -1;
    double filled_              = 0;
};

static_assert(ExecutionAlgo<VWAPAlgo>);

// ---------------------------------------------------------------------------
// IS — Implementation Shortfall (Almgren-Chriss optimal execution).
//
// Front-loads execution to reduce exposure to price drift at the cost of
// higher market impact.  The urgency parameter λ controls the trade-off:
//   λ → 0: trade everything immediately (max impact, zero drift risk)
//   λ → ∞: spread evenly (TWAP-like, zero impact, max drift risk)
//
// The optimal schedule is: q(t) = Q · sinh(κ(T-t)) / sinh(κT)
// where κ = √(λ/η) is the decay rate.
// ---------------------------------------------------------------------------
class ISAlgo {
public:
    explicit ISAlgo(ExecutionOrder order, int n_slices = 20,
                    double urgency = 1.0, double eta = 0.2)
        : order_(std::move(order)), n_slices_(n_slices)
    {
        // Build Almgren-Chriss optimal schedule
        double kappa = std::sqrt(urgency / std::max(eta, 1e-6));
        double T = 1.0;  // normalised horizon
        schedule_.resize(n_slices);
        double total = 0;
        for (int i = 0; i < n_slices; ++i) {
            double t = static_cast<double>(i) / n_slices;
            schedule_[i] = std::sinh(kappa * (T - t)) / std::sinh(kappa * T);
            total += schedule_[i];
        }
        // Normalise so total = order quantity
        for (auto& s : schedule_) s = (s / total) * order_.total_qty;

        slice_interval_us_ = (order_.end_time_us - order_.start_time_us) / n_slices;
        next_slice_time_ = order_.start_time_us + slice_interval_us_;
    }

    ISAlgo(ISAlgo&&) noexcept = default;
    ISAlgo& operator=(ISAlgo&&) noexcept = default;

    [[nodiscard]] std::optional<SliceDecision> on_tick(const MarketTick& tick) {
        if (filled_ >= order_.total_qty) return std::nullopt;
        if (tick.timestamp_us < next_slice_time_) return std::nullopt;

        int idx = std::min(slice_idx_, n_slices_ - 1);
        double qty = std::min(schedule_[idx], order_.total_qty - filled_);
        ++slice_idx_;
        next_slice_time_ += slice_interval_us_;
        filled_ += qty;

        return SliceDecision{qty, filled_ >= order_.total_qty};
    }

    [[nodiscard]] double filled() const noexcept { return filled_; }
    [[nodiscard]] std::string name() const { return "IS-Almgren"; }

private:
    ExecutionOrder order_;
    int    n_slices_;
    std::vector<double> schedule_;
    double slice_interval_us_ = 0;
    double next_slice_time_   = 0;
    int    slice_idx_         = 0;
    double filled_            = 0;
};

static_assert(ExecutionAlgo<ISAlgo>);

} // namespace pricer::execution
