#pragma once
/// @file performance.hpp
/// Performance analytics — Sharpe, Sortino, Calmar, max drawdown from any
/// PnL or returns series.
///
/// Why: After running a backtest you need portfolio-level risk/return metrics
/// to compare strategies, tune parameters, and report to stakeholders.
/// This module takes a raw PnL or cumulative-equity series and computes
/// annualized return, volatility, Sharpe, Sortino, Calmar, max drawdown
/// (value and duration).

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <string>
#include <format>
#include <concepts>
#include <ranges>
#include <limits>

namespace pricer::execution {

// ---------------------------------------------------------------------------
// ReturnSeries concept — any range whose elements convert to double.
// ---------------------------------------------------------------------------
template <typename R>
concept ReturnSeries = std::ranges::range<R> &&
    std::convertible_to<std::ranges::range_value_t<R>, double>;

// ---------------------------------------------------------------------------
// PerformanceConfig — parameters for annualization and benchmarking.
// ---------------------------------------------------------------------------
struct PerformanceConfig {
    double annualization_factor = 252.0;  // trading days per year
    double risk_free_rate       = 0.0;    // annualized risk-free rate
};

// ---------------------------------------------------------------------------
// PerformanceReport — all computed metrics in one struct.
// ---------------------------------------------------------------------------
struct PerformanceReport {
    double annualized_return  = 0.0;
    double annualized_vol     = 0.0;
    double sharpe_ratio       = 0.0;
    double max_drawdown       = 0.0;   // as positive fraction, e.g. 0.15 = 15%
    int    max_dd_duration    = 0;     // periods in longest drawdown
    double sortino_ratio      = 0.0;
    double calmar_ratio       = 0.0;
    int    n_periods          = 0;     // number of return observations
};

// ---------------------------------------------------------------------------
// cumulative_to_returns — convert a cumulative PnL / equity series to
// period-over-period returns.  Requires at least 2 data points.
// ---------------------------------------------------------------------------
inline std::vector<double> cumulative_to_returns(const std::vector<double>& cumulative) {
    if (cumulative.size() < 2) return {};
    std::vector<double> ret;
    ret.reserve(cumulative.size() - 1);
    for (std::size_t i = 1; i < cumulative.size(); ++i) {
        double prev = cumulative[i - 1];
        if (prev == 0.0) {
            ret.push_back(0.0);
        } else {
            ret.push_back((cumulative[i] - prev) / std::abs(prev));
        }
    }
    return ret;
}

// ---------------------------------------------------------------------------
// compute_performance — takes a period-returns series and config, returns
// a full PerformanceReport.
// ---------------------------------------------------------------------------
template <ReturnSeries R>
[[nodiscard]] PerformanceReport compute_performance(
    const R& returns,
    PerformanceConfig cfg = {}
) {
    PerformanceReport rpt;

    // Materialise into a vector for multi-pass
    std::vector<double> rets(std::ranges::begin(returns), std::ranges::end(returns));
    rpt.n_periods = static_cast<int>(rets.size());

    if (rets.empty()) return rpt;
    if (rets.size() == 1) {
        rpt.annualized_return = rets[0] * cfg.annualization_factor;
        rpt.annualized_vol    = 0.0;
        rpt.sharpe_ratio      = 0.0;
        rpt.sortino_ratio     = 0.0;
        return rpt;
    }

    const auto n = static_cast<double>(rets.size());

    // --- Mean period return ---
    double sum = std::accumulate(rets.begin(), rets.end(), 0.0);
    double mean_ret = sum / n;

    // --- Annualized return (geometric) ---
    // Compound: (1+r1)(1+r2)...(1+rn)^(AF/n) - 1
    double compound = 1.0;
    for (double r : rets) compound *= (1.0 + r);

    if (compound > 0.0) {
        rpt.annualized_return = std::pow(compound, cfg.annualization_factor / n) - 1.0;
    } else {
        // Total loss — annualized return is -100%
        rpt.annualized_return = -1.0;
    }

    // --- Annualized volatility ---
    double var_sum = 0.0;
    for (double r : rets) {
        double diff = r - mean_ret;
        var_sum += diff * diff;
    }
    double period_vol = std::sqrt(var_sum / (n - 1.0));
    rpt.annualized_vol = period_vol * std::sqrt(cfg.annualization_factor);

    // --- Sharpe ratio ---
    double rf_per_period = cfg.risk_free_rate / cfg.annualization_factor;
    double excess_mean = mean_ret - rf_per_period;
    rpt.sharpe_ratio = (period_vol > 1e-15)
        ? (excess_mean / period_vol) * std::sqrt(cfg.annualization_factor)
        : 0.0;

    // --- Sortino ratio (downside deviation) ---
    double downside_sum = 0.0;
    int    downside_count = 0;
    for (double r : rets) {
        double excess = r - rf_per_period;
        if (excess < 0.0) {
            downside_sum += excess * excess;
            ++downside_count;
        }
    }
    if (downside_count > 0) {
        double downside_dev = std::sqrt(downside_sum / n);
        rpt.sortino_ratio = (downside_dev > 1e-15)
            ? (excess_mean / downside_dev) * std::sqrt(cfg.annualization_factor)
            : std::numeric_limits<double>::infinity();
    } else {
        // No negative returns — Sortino is infinity (best case)
        rpt.sortino_ratio = std::numeric_limits<double>::infinity();
    }

    // --- Max drawdown (value and duration) ---
    double peak = 1.0;
    double equity = 1.0;
    double max_dd = 0.0;
    int    dd_start = 0;
    int    longest_dd = 0;
    bool   in_dd = false;

    for (int i = 0; i < static_cast<int>(rets.size()); ++i) {
        equity *= (1.0 + rets[static_cast<std::size_t>(i)]);
        if (equity > peak) {
            peak = equity;
            if (in_dd) {
                longest_dd = std::max(longest_dd, i - dd_start);
                in_dd = false;
            }
        } else {
            double dd = (peak - equity) / peak;
            if (dd > max_dd) max_dd = dd;
            if (!in_dd) {
                dd_start = i;
                in_dd = true;
            }
        }
    }
    // Close any open drawdown at the end
    if (in_dd) {
        longest_dd = std::max(longest_dd, static_cast<int>(rets.size()) - dd_start);
    }

    rpt.max_drawdown    = max_dd;
    rpt.max_dd_duration = longest_dd;

    // --- Calmar ratio ---
    rpt.calmar_ratio = (max_dd > 1e-15)
        ? rpt.annualized_return / max_dd
        : 0.0;

    return rpt;
}

// ---------------------------------------------------------------------------
// format_report — human-readable summary of a PerformanceReport.
// ---------------------------------------------------------------------------
[[nodiscard]] inline std::string format_report(const PerformanceReport& r) {
    auto fmt_sortino = [&]() -> std::string {
        if (std::isinf(r.sortino_ratio))
            return "       +Inf";
        return std::format("{:>+11.4f}", r.sortino_ratio);
    };

    return std::format(
        "Performance Report ({} periods)\n"
        "  Annualized return:  {:>+11.4f}  ({:+.2f}%)\n"
        "  Annualized vol:     {:>11.4f}  ({:.2f}%)\n"
        "  Sharpe ratio:       {:>+11.4f}\n"
        "  Sortino ratio:      {}\n"
        "  Max drawdown:       {:>11.4f}  ({:.2f}%)\n"
        "  Max DD duration:    {:>11d}  periods\n"
        "  Calmar ratio:       {:>+11.4f}",
        r.n_periods,
        r.annualized_return, r.annualized_return * 100,
        r.annualized_vol, r.annualized_vol * 100,
        r.sharpe_ratio,
        fmt_sortino(),
        r.max_drawdown, r.max_drawdown * 100,
        r.max_dd_duration,
        r.calmar_ratio);
}

} // namespace pricer::execution
