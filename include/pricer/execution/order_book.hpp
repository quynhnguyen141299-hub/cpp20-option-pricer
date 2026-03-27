#pragma once
/// @file order_book.hpp
/// Simulated limit order book for backtesting execution algorithms.
///
/// Why: Algorithmic execution quality depends on market microstructure —
/// bid-ask spread, queue position, fill probability, and market impact.
/// This LOB simulator provides a lightweight but realistic environment
/// to test TWAP/VWAP/IS algos before deploying to a real venue.
///
/// Design: Level-2 book with configurable depth, spread dynamics, and
/// a square-root market impact model (empirically validated by
/// Almgren-Chriss and JPM's own research).

#include "../core/types.hpp"
#include "../core/errors.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <format>
#include <chrono>
#include <cstdint>
#include <string>

namespace pricer::execution {

// ---------------------------------------------------------------------------
// Side — bid or ask.
// ---------------------------------------------------------------------------
enum class Side : std::uint8_t { Bid, Ask };

[[nodiscard]] constexpr std::string_view str(Side s) noexcept {
    return s == Side::Bid ? "Bid" : "Ask";
}

// ---------------------------------------------------------------------------
// PriceLevel — one level in the order book.
// ---------------------------------------------------------------------------
struct PriceLevel {
    double price    = 0;
    double quantity = 0;   // in base currency units (e.g. EUR for EURUSD)
};

// ---------------------------------------------------------------------------
// MarketTick — a point-in-time snapshot of the top of book.
// ---------------------------------------------------------------------------
struct MarketTick {
    double timestamp_us = 0;  // microseconds since epoch (or sim start)
    double mid          = 0;
    double bid          = 0;
    double ask          = 0;
    double spread       = 0;
    double bid_qty      = 0;
    double ask_qty      = 0;
    double last_trade   = 0;
    double volume       = 0;  // cumulative session volume

    [[nodiscard]] double spread_bps() const noexcept {
        return mid > 0 ? (spread / mid) * 10'000.0 : 0;
    }
};

// ---------------------------------------------------------------------------
// Fill — the result of executing a slice against the book.
// ---------------------------------------------------------------------------
struct Fill {
    double price        = 0;
    double quantity     = 0;
    double slippage_bps = 0;  // vs mid at time of execution
    double timestamp_us = 0;
    Side   side         = Side::Bid;

    [[nodiscard]] std::string str() const {
        return std::format("{} {:.0f} @ {:.5f} (slip={:+.2f}bps)",
                           execution::str(side), quantity, price, slippage_bps);
    }
};

// ---------------------------------------------------------------------------
// ImpactModel — square-root temporary market impact.
//
// ΔP/P = η · σ · √(q / V_daily)
//
// where η is the impact coefficient (typically 0.1–0.5 for FX),
// σ is daily vol, q is order size, V_daily is daily volume.
// This is the Almgren-Chriss model used on most execution desks.
// ---------------------------------------------------------------------------
struct ImpactModel {
    double eta          = 0.2;     // impact coefficient
    double daily_vol    = 0.007;   // daily σ (≈ 7.5% annualised / √252)
    double daily_volume = 50e9;    // typical EURUSD daily volume

    /// Temporary impact in price units for a given order size.
    [[nodiscard]] double impact(double order_size, double mid) const noexcept {
        double participation = std::abs(order_size) / daily_volume;
        return eta * daily_vol * mid * std::sqrt(participation);
    }

    /// Impact in basis points.
    [[nodiscard]] double impact_bps(double order_size, double mid) const noexcept {
        return (impact(order_size, mid) / mid) * 10'000.0;
    }
};

// ---------------------------------------------------------------------------
// OrderBookSim — generates realistic tick data and simulates fills.
//
// RAII: owns its RNG state.  Move-only (no accidental RNG correlation).
// ---------------------------------------------------------------------------
struct OrderBookConfig {
    double initial_mid    = 1.0850;
    double half_spread    = 0.00005;  // 0.5 pips for EURUSD
    double tick_size      = 0.00001;  // 1 pip = 0.0001, tick = 0.1 pip
    double vol_per_tick   = 0.000002; // mid-price random walk step
    double base_qty       = 5e6;      // base quantity at each level
    int    book_depth     = 5;        // levels per side
    double daily_volume   = 50e9;
};

class OrderBookSim {
public:
    using Config = OrderBookConfig;

    explicit OrderBookSim(Config cfg = {}, std::uint64_t seed = 42)
        : cfg_(cfg), rng_(seed), mid_(cfg.initial_mid), cumulative_vol_(0)
    {
        rebuild_book();
    }

    // Move-only
    OrderBookSim(OrderBookSim&&) noexcept = default;
    OrderBookSim& operator=(OrderBookSim&&) noexcept = default;
    OrderBookSim(const OrderBookSim&) = delete;
    OrderBookSim& operator=(const OrderBookSim&) = delete;

    /// Advance the book by one tick (random walk + spread noise).
    MarketTick step(double timestamp_us) noexcept {
        // Random walk on mid
        std::normal_distribution<double> mid_dist(0.0, cfg_.vol_per_tick);
        mid_ += mid_dist(rng_);

        // Slight spread variation
        std::uniform_real_distribution<double> spread_noise(0.8, 1.2);
        double half = cfg_.half_spread * spread_noise(rng_);

        // Volume increment
        std::exponential_distribution<double> vol_dist(1.0 / 1e6);
        cumulative_vol_ += vol_dist(rng_);

        rebuild_book();

        return MarketTick{
            .timestamp_us = timestamp_us,
            .mid          = mid_,
            .bid          = mid_ - half,
            .ask          = mid_ + half,
            .spread       = 2.0 * half,
            .bid_qty      = bids_.empty() ? 0 : bids_[0].quantity,
            .ask_qty      = asks_.empty() ? 0 : asks_[0].quantity,
            .last_trade   = mid_,
            .volume       = cumulative_vol_
        };
    }

    /// Simulate a fill against the book with market impact.
    [[nodiscard]] Fill execute(double quantity, Side side,
                               const ImpactModel& impact_model,
                               double timestamp_us) noexcept {
        double raw_price = (side == Side::Ask) ? asks_[0].price : bids_[0].price;
        double imp = impact_model.impact(quantity, mid_);

        // Buy → price moves up; sell → price moves down
        double fill_price = (side == Side::Ask)
            ? raw_price + imp
            : raw_price - imp;

        double slip = (fill_price - mid_) / mid_ * 10'000.0;
        if (side == Side::Bid) slip = -slip;  // selling: slip is negative if worse

        // Adjust mid post-impact (temporary impact decays, but shift mid slightly)
        mid_ += (side == Side::Ask ? 1 : -1) * imp * 0.3;  // 30% permanent

        return Fill{
            .price        = fill_price,
            .quantity     = quantity,
            .slippage_bps = slip,
            .timestamp_us = timestamp_us,
            .side         = side
        };
    }

    [[nodiscard]] double mid() const noexcept { return mid_; }
    [[nodiscard]] const std::vector<PriceLevel>& bids() const noexcept { return bids_; }
    [[nodiscard]] const std::vector<PriceLevel>& asks() const noexcept { return asks_; }
    [[nodiscard]] const Config& config() const noexcept { return cfg_; }

private:
    void rebuild_book() {
        bids_.resize(cfg_.book_depth);
        asks_.resize(cfg_.book_depth);

        std::uniform_real_distribution<double> qty_noise(0.5, 1.5);
        for (int i = 0; i < cfg_.book_depth; ++i) {
            bids_[i] = {mid_ - cfg_.half_spread - i * cfg_.tick_size * 10,
                         cfg_.base_qty * qty_noise(rng_)};
            asks_[i] = {mid_ + cfg_.half_spread + i * cfg_.tick_size * 10,
                         cfg_.base_qty * qty_noise(rng_)};
        }
    }

    Config cfg_;
    std::mt19937_64 rng_;
    double mid_;
    double cumulative_vol_;
    std::vector<PriceLevel> bids_;
    std::vector<PriceLevel> asks_;
};

} // namespace pricer::execution
