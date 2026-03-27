#pragma once
/// @file tca.hpp
/// Transaction Cost Analysis (TCA) — post-trade execution quality metrics.
///
/// Why: Every execution desk measures algo performance after the fact.
/// The key benchmarks are:
///   - Arrival price shortfall: how much worse is the avg fill vs the mid
///     at the moment the order arrived?  (The IS benchmark.)
///   - VWAP slippage: how does the avg fill compare to the session VWAP?
///   - Spread capture: did the algo fill inside the spread or cross it?
///   - Market impact: what was the total cost attributable to our order?
///
/// These metrics feed back into algo parameter tuning (urgency, slice size).

#include "order_book.hpp"
#include <vector>
#include <numeric>
#include <cmath>
#include <format>
#include <string>

namespace pricer::execution {

// ---------------------------------------------------------------------------
// TCAReport — comprehensive execution quality analysis.
// ---------------------------------------------------------------------------
struct TCAReport {
    std::string algo_name;

    // Order summary
    double total_qty       = 0;
    int    n_fills         = 0;
    double duration_us     = 0;  // first fill to last fill

    // Price benchmarks
    double arrival_mid     = 0;  // mid at order arrival
    double avg_fill_price  = 0;  // volume-weighted avg fill
    double vwap_benchmark  = 0;  // market VWAP over execution window
    double terminal_mid    = 0;  // mid at end of execution

    // Cost metrics (in basis points)
    double arrival_cost_bps  = 0;  // vs arrival mid (IS benchmark)
    double vwap_slip_bps     = 0;  // vs market VWAP
    double spread_cost_bps   = 0;  // avg half-spread paid
    double impact_cost_bps   = 0;  // total market impact

    // Risk metrics
    double timing_risk_bps   = 0;  // mid drift from arrival to terminal
    double fill_rate         = 0;  // fraction of order filled
    double participation_rate = 0; // our volume / market volume

    [[nodiscard]] std::string summary() const {
        return std::format(
            "[{}] TCA Report\n"
            "  Fills: {} | Qty: {:.0f} | Duration: {:.0f}ms | Fill rate: {:.1f}%\n"
            "  Avg fill: {:.5f} | Arrival mid: {:.5f} | Terminal mid: {:.5f}\n"
            "  ── Cost Breakdown (bps) ──\n"
            "  Arrival shortfall:  {:+.2f}\n"
            "  VWAP slippage:      {:+.2f}\n"
            "  Spread cost:        {:+.2f}\n"
            "  Impact cost:        {:+.2f}\n"
            "  ── Risk ──\n"
            "  Timing risk:        {:+.2f} bps\n"
            "  Participation:      {:.1f}%",
            algo_name, n_fills, total_qty, duration_us / 1000.0, fill_rate * 100,
            avg_fill_price, arrival_mid, terminal_mid,
            arrival_cost_bps, vwap_slip_bps, spread_cost_bps, impact_cost_bps,
            timing_risk_bps, participation_rate * 100);
    }
};

// ---------------------------------------------------------------------------
// TCAEngine — computes TCA from a sequence of fills and market ticks.
// ---------------------------------------------------------------------------
class TCAEngine {
public:
    /// Compute TCA from fills and the tick stream that was active during execution.
    [[nodiscard]] static TCAReport analyse(
        const std::string& algo_name,
        const std::vector<Fill>& fills,
        const std::vector<MarketTick>& ticks,
        double arrival_mid,
        double total_order_qty,
        Side side
    ) {
        TCAReport r;
        r.algo_name   = algo_name;
        r.arrival_mid = arrival_mid;
        r.total_qty   = total_order_qty;
        r.n_fills     = static_cast<int>(fills.size());

        if (fills.empty() || ticks.empty()) return r;

        // Volume-weighted average fill price
        double vwap_num = 0, vwap_den = 0;
        for (const auto& f : fills) {
            vwap_num += f.price * f.quantity;
            vwap_den += f.quantity;
        }
        r.avg_fill_price = vwap_den > 0 ? vwap_num / vwap_den : 0;

        // Market VWAP (from tick stream — proxy using mid * volume increments)
        double mkt_vwap_num = 0, mkt_vwap_den = 0;
        double prev_vol = 0;
        for (const auto& t : ticks) {
            double vol_inc = t.volume - prev_vol;
            if (vol_inc > 0) {
                mkt_vwap_num += t.mid * vol_inc;
                mkt_vwap_den += vol_inc;
            }
            prev_vol = t.volume;
        }
        r.vwap_benchmark = mkt_vwap_den > 0 ? mkt_vwap_num / mkt_vwap_den : arrival_mid;

        // Terminal mid
        r.terminal_mid = ticks.back().mid;

        // Duration
        r.duration_us = fills.back().timestamp_us - fills.front().timestamp_us;

        // Fill rate
        double filled = std::accumulate(fills.begin(), fills.end(), 0.0,
            [](double s, const Fill& f) { return s + f.quantity; });
        r.fill_rate = filled / total_order_qty;

        // Participation rate
        double market_vol = ticks.back().volume - ticks.front().volume;
        r.participation_rate = market_vol > 0 ? filled / market_vol : 0;

        // Cost metrics — sign convention: positive = cost (worse than benchmark)
        double sign = (side == Side::Ask) ? 1.0 : -1.0;  // buying: higher is worse

        // Arrival shortfall
        r.arrival_cost_bps = sign * (r.avg_fill_price - arrival_mid) / arrival_mid * 10'000.0;

        // VWAP slippage
        r.vwap_slip_bps = sign * (r.avg_fill_price - r.vwap_benchmark) / r.vwap_benchmark * 10'000.0;

        // Spread cost (average slippage from fills)
        double avg_slip = 0;
        for (const auto& f : fills) avg_slip += std::abs(f.slippage_bps);
        r.spread_cost_bps = fills.empty() ? 0 : avg_slip / fills.size();

        // Impact cost — difference between our avg fill and what a zero-impact
        // algo would have achieved (approximated by market VWAP)
        r.impact_cost_bps = sign * (r.avg_fill_price - r.vwap_benchmark) / r.vwap_benchmark * 10'000.0;

        // Timing risk — mid drift over execution window
        r.timing_risk_bps = sign * (r.terminal_mid - arrival_mid) / arrival_mid * 10'000.0;

        return r;
    }
};

} // namespace pricer::execution
