#pragma once
/// @file backtest.hpp
/// Backtest harness — drives an execution algo over simulated tick data
/// and produces a TCA report.
///
/// This is the glue between the order book simulator, the execution algo,
/// and the TCA engine.  It simulates the full lifecycle:
///   1. Generate tick stream from OrderBookSim
///   2. Feed each tick to the algo
///   3. When the algo wants to trade, execute against the book
///   4. Collect all fills
///   5. Run TCA on the fills + tick stream

#include "order_book.hpp"
#include "algo.hpp"
#include "tca.hpp"
#include <vector>
#include <format>
#include <iostream>

namespace pricer::execution {

// ---------------------------------------------------------------------------
// BacktestConfig
// ---------------------------------------------------------------------------
struct BacktestConfig {
    int    n_ticks          = 10'000;   // total ticks to simulate
    double tick_interval_us = 100'000;  // 100ms between ticks (10 ticks/sec)
    bool   verbose          = false;
};

// ---------------------------------------------------------------------------
// BacktestResult
// ---------------------------------------------------------------------------
struct BacktestResult {
    TCAReport               tca;
    std::vector<Fill>       fills;
    std::vector<MarketTick> ticks;
};

// ---------------------------------------------------------------------------
// run_backtest — generic over any ExecutionAlgo.
//
// The algo concept constraint means this function works with TWAP, VWAP,
// IS, or any future algo that satisfies the interface — without virtual
// dispatch, without inheritance hierarchies.
// ---------------------------------------------------------------------------
template <ExecutionAlgo Algo>
[[nodiscard]] BacktestResult run_backtest(
    Algo& algo,
    ExecutionOrder order,
    OrderBookSim& book,
    ImpactModel impact = {},
    BacktestConfig cfg = {}
) {
    std::vector<MarketTick> ticks;
    std::vector<Fill> fills;
    ticks.reserve(cfg.n_ticks);

    double arrival_mid = book.mid();

    for (int i = 0; i < cfg.n_ticks; ++i) {
        double ts = order.start_time_us + i * cfg.tick_interval_us;
        auto tick = book.step(ts);
        ticks.push_back(tick);

        auto decision = algo.on_tick(tick);
        if (decision && decision->quantity > 0) {
            auto fill = book.execute(decision->quantity, order.side, impact, ts);
            fills.push_back(fill);

            if (cfg.verbose) {
                std::cout << std::format("  t={:.0f}ms  {}\n",
                    ts / 1000.0, fill.str());
            }

            if (decision->is_final) break;
        }
    }

    auto tca = TCAEngine::analyse(
        algo.name(), fills, ticks, arrival_mid, order.total_qty, order.side);

    return BacktestResult{
        .tca   = std::move(tca),
        .fills = std::move(fills),
        .ticks = std::move(ticks)
    };
}

} // namespace pricer::execution
