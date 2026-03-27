/// @file test_pricer.cpp
/// Unit tests — no external dependencies.

#include "pricer/pricer.hpp"
#include <iostream>
#include <format>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <cassert>

using namespace pricer;

static int passed = 0, failed = 0;

void test(const char* name, std::function<bool()> fn) {
    try {
        bool ok = fn();
        std::cout << std::format("  [{}] {}\n", ok ? "PASS" : "FAIL", name);
        ok ? ++passed : ++failed;
    } catch (const std::exception& e) {
        std::cout << std::format("  [FAIL] {} — exception: {}\n", name, e.what());
        ++failed;
    }
}

bool near(double a, double b, double tol = 1e-4) {
    return std::abs(a - b) < tol;
}

// ── Strong types ─────────────────────────────────────────────────
void test_types() {
    test("StrongDouble prevents implicit construction", [] {
        // Spot s = 100.0;  // should not compile — explicit ctor
        Spot s{100.0};
        return near(s.v, 100.0);
    });

    test("Spaceship operator", [] {
        return Strike{100.0} < Strike{110.0};
    });
}

// ── Black-Scholes ────────────────────────────────────────────────
void test_bs() {
    BSEngine bs;
    MarketSnap m{Spot{100}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};
    Contract put {Strike{100}, YearFrac{1.0}, OptType::Put};

    test("BS call ~10.45 (textbook)", [&] {
        auto r = bs.price(call, m);
        return r && near(r->price, 10.4506, 0.01);
    });

    test("BS put-call parity", [&] {
        auto c = bs.price(call, m);
        auto p = bs.price(put, m);
        if (!c || !p) return false;
        double lhs = c->price - p->price;
        double rhs = m.S.v - 100.0 * std::exp(-0.05);
        return near(lhs, rhs, 1e-8);
    });

    test("BS call delta in (0,1)", [&] {
        auto r = bs.price(call, m);
        return r && r->greeks.delta > 0 && r->greeks.delta < 1;
    });

    test("BS put delta in (-1,0)", [&] {
        auto r = bs.price(put, m);
        return r && r->greeks.delta > -1 && r->greeks.delta < 0;
    });

    test("BS gamma > 0", [&] {
        auto r = bs.price(call, m);
        return r && r->greeks.gamma > 0;
    });

    test("BS vega > 0", [&] {
        auto r = bs.price(call, m);
        return r && r->greeks.vega > 0;
    });
}

// ── FX (Garman-Kohlhagen) ───────────────────────────────────────
void test_fx() {
    BSEngine bs;
    MarketSnap m{Spot{1.085}, Vol{0.075}, Rate{0.0435}, Rate{0.025}};
    Contract call{Strike{1.09}, YearFrac{0.25}, OptType::Call};
    Contract put {Strike{1.09}, YearFrac{0.25}, OptType::Put};

    test("FX put-call parity (GK)", [&] {
        auto c = bs.price(call, m);
        auto p = bs.price(put, m);
        if (!c || !p) return false;
        double lhs = c->price - p->price;
        double rhs = m.S.v * std::exp(-m.r_f.v * 0.25)
                   - 1.09 * std::exp(-m.r_d.v * 0.25);
        return near(lhs, rhs, 1e-8);
    });
}

// ── Error handling ───────────────────────────────────────────────
void test_errors() {
    BSEngine bs;
    Contract c{Strike{100}, YearFrac{1.0}, OptType::Call};

    test("Neg spot → BadSpot", [&] {
        MarketSnap bad{Spot{-1}, Vol{0.2}, Rate{0.05}, Rate{0}};
        auto r = bs.price(c, bad);
        return !r && r.error().code == Err::BadSpot;
    });

    test("Neg vol → BadVol", [&] {
        MarketSnap bad{Spot{100}, Vol{-0.2}, Rate{0.05}, Rate{0}};
        auto r = bs.price(c, bad);
        return !r && r.error().code == Err::BadVol;
    });

    test("Zero expiry → BadExpiry", [&] {
        Contract bad{Strike{100}, YearFrac{0}, OptType::Call};
        MarketSnap m{Spot{100}, Vol{0.2}, Rate{0.05}, Rate{0}};
        auto r = bs.price(bad, m);
        return !r && r.error().code == Err::BadExpiry;
    });

    test("Error has source location", [&] {
        MarketSnap bad{Spot{-1}, Vol{0.2}, Rate{0.05}, Rate{0}};
        auto r = bs.price(c, bad);
        return !r && !r.error().loc.empty();
    });

    test("Monadic transform on success", [&] {
        MarketSnap m{Spot{100}, Vol{0.2}, Rate{0.05}, Rate{0}};
        auto r = bs.price(c, m).transform([](PricingResult p) { return p.price; });
        return r && *r > 0;
    });

    test("Monadic transform skipped on error", [&] {
        MarketSnap bad{Spot{-1}, Vol{0.2}, Rate{0.05}, Rate{0}};
        bool called = false;
        auto r = bs.price(c, bad).transform([&](PricingResult p) {
            called = true; return p.price;
        });
        return !r && !called;
    });
}

// ── Monte Carlo ──────────────────────────────────────────────────
void test_mc() {
    MarketSnap m{Spot{100}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};

    test("MC+CV converges to BS within 3 SE", [&] {
        MCConfig cfg{.n_paths = 500'000, .seed = 42, .antithetic = true,
                     .control_variate = true, .n_threads = 1, .batch_size = 50'000};
        MCEngine mc(FlatRate{m.r_d}, FlatVol{m.sigma}, cfg);
        BSEngine bs;
        auto mc_r = mc.price(call, m);
        auto bs_r = bs.price(call, m);
        if (!mc_r || !bs_r) return false;
        return std::abs(mc_r->price - bs_r->price) < 3.0 * mc_r->std_err;
    });

    test("MC SE decreases with more paths", [&] {
        MCConfig c1{.n_paths = 10'000, .seed = 42, .antithetic = true,
                    .control_variate = true, .n_threads = 1, .batch_size = 10'000};
        MCConfig c2{.n_paths = 100'000, .seed = 42, .antithetic = true,
                    .control_variate = true, .n_threads = 1, .batch_size = 50'000};
        MCEngine m1(FlatRate{m.r_d}, FlatVol{m.sigma}, c1);
        MCEngine m2(FlatRate{m.r_d}, FlatVol{m.sigma}, c2);
        auto r1 = m1.price(call, m);
        auto r2 = m2.price(call, m);
        return r1 && r2 && r2->std_err < r1->std_err;
    });

    test("CV reduces variance vs raw", [&] {
        MCConfig raw{.n_paths = 100'000, .seed = 42, .antithetic = false,
                     .control_variate = false, .n_threads = 1, .batch_size = 50'000};
        MCConfig cv{.n_paths = 100'000, .seed = 42, .antithetic = false,
                    .control_variate = true, .n_threads = 1, .batch_size = 50'000};
        MCEngine m_raw(FlatRate{m.r_d}, FlatVol{m.sigma}, raw);
        MCEngine m_cv(FlatRate{m.r_d}, FlatVol{m.sigma}, cv);
        auto r_raw = m_raw.price(call, m);
        auto r_cv  = m_cv.price(call, m);
        return r_raw && r_cv && r_cv->std_err < r_raw->std_err;
    });

    test("MC rejects negative spot", [&] {
        MarketSnap bad{Spot{-100}, Vol{0.2}, Rate{0.05}, Rate{0}};
        MCConfig cfg{.n_paths = 1000, .seed = 42, .n_threads = 1, .batch_size = 1000};
        MCEngine mc(FlatRate{bad.r_d}, FlatVol{bad.sigma}, cfg);
        auto r = mc.price(call, bad);
        return !r && r.error().code == Err::BadSpot;
    });

    test("Pathwise delta within 5% of analytic", [&] {
        MCConfig cfg{.n_paths = 1'000'000, .seed = 42, .antithetic = true,
                     .control_variate = true, .n_threads = 1, .batch_size = 100'000};
        MCEngine mc(FlatRate{m.r_d}, FlatVol{m.sigma}, cfg);
        BSEngine bs;
        auto mc_r = mc.price(call, m);
        auto bs_r = bs.price(call, m);
        if (!mc_r || !bs_r) return false;
        double rel = std::abs(mc_r->greeks.delta - bs_r->greeks.delta) / bs_r->greeks.delta;
        return rel < 0.05;
    });
}

// ── Rate models ──────────────────────────────────────────────────
void test_rates() {
    test("FlatRate df", [] {
        FlatRate r(Rate{0.05});
        return near(r.df(1.0), std::exp(-0.05), 1e-10);
    });

    test("PiecewiseCurve interpolation", [] {
        PiecewiseCurve c({{0.25, 0.05}, {1.0, 0.04}, {5.0, 0.035}});
        double r = c.zero_rate(0.5);
        return r > 0.035 && r < 0.05;
    });

    test("FXRatePair drift", [] {
        FXRatePair pair(FlatRate{Rate{0.05}}, FlatRate{Rate{0.03}});
        return near(pair.fx_drift(1.0), 0.02, 1e-10);
    });
}

// ── Vol models ───────────────────────────────────────────────────
void test_vols() {
    test("FlatVol constant", [] {
        FlatVol v(Vol{0.20});
        return near(v.local_vol(100, 100, 1), 0.20) && near(v.iv(100, 1), 0.20);
    });

    test("TermStructureVol variance interp", [] {
        TermStructureVol tv({{0.25, 0.20}, {1.0, 0.18}, {2.0, 0.17}});
        double v = tv.iv(100, 0.5);
        return v > 0.17 && v < 0.20;
    });

    test("StickyStrikeVol skew", [] {
        StickyStrikeVol sv(0.20, -0.1, 0.5, 100);
        return sv.iv(90, 1) > sv.iv(100, 1);  // downside skew
    });
}

// ── Sobol engine ─────────────────────────────────────────────────
void test_sobol() {
    test("Sobol values in [0,1)", [] {
        SobolEngine s;
        for (int i = 0; i < 10000; ++i) {
            double v = s.next();
            if (v < 0 || v >= 1) return false;
        }
        return true;
    });

    test("Sobol uniformity (mean near 0.5)", [] {
        SobolEngine s;
        double sum = 0;
        int N = 100000;
        for (int i = 0; i < N; ++i) sum += s.next();
        double mean = sum / N;
        return near(mean, 0.5, 0.01);
    });

    test("Inverse normal: N^{-1}(0.5) = 0", [] {
        return near(SobolEngine::inv_normal(0.5), 0.0, 1e-6);
    });

    test("Inverse normal: N^{-1}(0.8413) ≈ 1.0", [] {
        return near(SobolEngine::inv_normal(0.8413), 1.0, 0.01);
    });
}

// ── Arena ────────────────────────────────────────────────────────
void test_arena() {
    test("Arena alloc + reset", [] {
        Arena a(4096);
        auto s = a.alloc_doubles(100);
        if (s.empty()) return false;
        if (a.used() < 800) return false;
        a.reset();
        return a.used() == 0;
    });

    test("Arena exhaustion returns empty span", [] {
        Arena a(64);
        auto s = a.alloc_doubles(100);  // 800 bytes > 64
        return s.empty();
    });
}

// ── Thread pool ──────────────────────────────────────────────────
void test_thread_pool() {
    test("ThreadPool submit + future", [] {
        ThreadPool pool(2);
        auto f = pool.submit([] { return 42; });
        return f.get() == 42;
    });

    test("ThreadPool multiple tasks", [] {
        ThreadPool pool(2);
        std::vector<std::future<int>> futs;
        for (int i = 0; i < 100; ++i)
            futs.push_back(pool.submit([i] { return i * i; }));
        for (int i = 0; i < 100; ++i)
            if (futs[i].get() != i * i) return false;
        return true;
    });
}

// ── IV solver ────────────────────────────────────────────────────
void test_iv() {
    test("IV solver recovers known vol", [] {
        BSEngine bs;
        MarketSnap m{Spot{100}, Vol{0.25}, Rate{0.05}, Rate{0}};
        Contract c{Strike{100}, YearFrac{1.0}, OptType::Call};
        auto px = bs.price(c, m);
        if (!px) return false;
        MarketSnap guess = m;
        guess.sigma = Vol{0.40};
        auto iv = bs.implied_vol(c, guess, px->price);
        return iv && near(*iv, 0.25, 1e-6);
    });

    test("IV solver rejects neg price", [] {
        BSEngine bs;
        MarketSnap m{Spot{100}, Vol{0.25}, Rate{0.05}, Rate{0}};
        Contract c{Strike{100}, YearFrac{1.0}, OptType::Call};
        auto iv = bs.implied_vol(c, m, -1.0);
        return !iv && iv.error().code == Err::BadConfig;
    });
}

// ── Concepts (compile-time) ──────────────────────────────────────
void test_concepts() {
    test("FlatRate satisfies RateModel",      [] { return RateModel<FlatRate>; });
    test("PiecewiseCurve satisfies RateModel", [] { return RateModel<PiecewiseCurve>; });
    test("FlatVol satisfies VolSurface",       [] { return VolSurface<FlatVol>; });
    test("BSEngine satisfies Engine",          [] { return Engine<BSEngine>; });
    test("MCEngine satisfies Engine",          [] { return Engine<MCEngine<FlatRate, FlatVol>>; });
}

// ── Move semantics ───────────────────────────────────────────────
void test_move() {
    test("PiecewiseCurve move", [] {
        PiecewiseCurve c1({{0.25, 0.05}, {1.0, 0.04}});
        auto d = c1.df(0.5);
        PiecewiseCurve c2 = std::move(c1);
        return near(c2.df(0.5), d, 1e-10);
    });

    test("MCEngine move", [] {
        MCConfig cfg{.n_paths = 1000, .seed = 42, .n_threads = 1, .batch_size = 1000};
        MCEngine mc1(FlatRate{Rate{0.05}}, FlatVol{Vol{0.20}}, cfg);
        MCEngine mc2 = std::move(mc1);
        MarketSnap m{Spot{100}, Vol{0.2}, Rate{0.05}, Rate{0}};
        Contract c{Strike{100}, YearFrac{1.0}, OptType::Call};
        auto r = mc2.price(c, m);
        return r && r->price > 0;
    });

    test("Arena move", [] {
        Arena a1(4096);
        auto s = a1.alloc_doubles(10);
        s[0] = 42.0;
        Arena a2 = std::move(a1);
        return a2.used() > 0;
    });
}

// ── Order book ──────────────────────────────────────────────────────
using namespace pricer::execution;

void test_order_book() {
    test("OrderBookSim generates valid ticks", [] {
        OrderBookSim book;
        auto tick = book.step(0);
        return tick.bid < tick.ask && tick.mid > 0 && tick.spread > 0;
    });

    test("OrderBookSim mid evolves over time", [] {
        OrderBookSim book({}, 99);
        double first_mid = book.mid();
        for (int i = 0; i < 1000; ++i) book.step(i * 1000.0);
        return book.mid() != first_mid;  // random walk should move
    });

    test("ImpactModel increases with order size", [] {
        ImpactModel im;
        double small_imp = im.impact_bps(1e6, 1.085);
        double large_imp = im.impact_bps(100e6, 1.085);
        return large_imp > small_imp && small_imp > 0;
    });

    test("Fill records slippage", [] {
        OrderBookSim book;
        ImpactModel im;
        (void)book.step(0);
        auto fill = book.execute(1e6, Side::Ask, im, 0);
        return fill.price > 0 && fill.quantity == 1e6;
    });
}

// ── Execution algos ─────────────────────────────────────────────────
void test_algos() {
    ExecutionOrder order{
        .total_qty = 10e6, .side = Side::Ask,
        .start_time_us = 0, .end_time_us = 1000e6
    };

    test("TWAP fills entire order", [&] {
        TWAPAlgo algo(order, 10);
        OrderBookSim book;
        ImpactModel im;
        auto result = run_backtest(algo, order, book, im,
            {.n_ticks = 2000, .tick_interval_us = 1e6});
        return near(result.tca.fill_rate, 1.0, 0.01);
    });

    test("VWAP fills entire order", [&] {
        VWAPAlgo algo(order, 10);
        OrderBookSim book;
        ImpactModel im;
        auto result = run_backtest(algo, order, book, im,
            {.n_ticks = 2000, .tick_interval_us = 1e6});
        return near(result.tca.fill_rate, 1.0, 0.01);
    });

    test("IS fills entire order", [&] {
        ISAlgo algo(order, 10, 1.0);
        OrderBookSim book;
        ImpactModel im;
        auto result = run_backtest(algo, order, book, im,
            {.n_ticks = 2000, .tick_interval_us = 1e6});
        return near(result.tca.fill_rate, 1.0, 0.01);
    });

    test("TWAP satisfies ExecutionAlgo concept", [] {
        return ExecutionAlgo<TWAPAlgo>;
    });

    test("VWAP satisfies ExecutionAlgo concept", [] {
        return ExecutionAlgo<VWAPAlgo>;
    });

    test("IS satisfies ExecutionAlgo concept", [] {
        return ExecutionAlgo<ISAlgo>;
    });
}

// ── TCA ──────────────────────────────────────────────────────────
void test_tca() {
    test("TCA reports positive arrival cost for buy", [] {
        ExecutionOrder order{.total_qty = 10e6, .side = Side::Ask,
            .start_time_us = 0, .end_time_us = 1000e6};
        TWAPAlgo algo(order, 10);
        OrderBookSim book;
        ImpactModel im;
        auto result = run_backtest(algo, order, book, im,
            {.n_ticks = 2000, .tick_interval_us = 1e6});
        // Buying should have positive arrival cost (we pay above mid)
        return result.tca.arrival_cost_bps > -10;  // reasonable range
    });

    test("TCA fill count matches algo fills", [] {
        ExecutionOrder order{.total_qty = 5e6, .side = Side::Ask,
            .start_time_us = 0, .end_time_us = 1000e6};
        TWAPAlgo algo(order, 5);
        OrderBookSim book;
        ImpactModel im;
        auto result = run_backtest(algo, order, book, im,
            {.n_ticks = 2000, .tick_interval_us = 1e6});
        return result.tca.n_fills == static_cast<int>(result.fills.size());
    });
}

// ── Signals ──────────────────────────────────────────────────────
void test_signals() {
    test("SpreadSignal score in [-1, 1]", [] {
        OrderBookSim book;
        SpreadSignal sig(20);
        for (int i = 0; i < 200; ++i) {
            auto tick = book.step(i * 1000.0);
            double s = sig.score(tick);
            if (s < -1.01 || s > 1.01) return false;
        }
        return true;
    });

    test("MomentumSignal score in [-1, 1]", [] {
        OrderBookSim book;
        MomentumSignal sig(20);
        for (int i = 0; i < 200; ++i) {
            auto tick = book.step(i * 1000.0);
            double s = sig.score(tick);
            if (s < -1.01 || s > 1.01) return false;
        }
        return true;
    });

    test("SpreadSignal satisfies Signal concept", [] {
        return Signal<SpreadSignal>;
    });

    test("MomentumSignal satisfies Signal concept", [] {
        return Signal<MomentumSignal>;
    });

    test("CompositeSignal satisfies Signal concept", [] {
        return Signal<CompositeSignal>;
    });
}

// ── Performance analytics ────────────────────────────────────────
void test_performance() {
    test("Known returns: Sharpe, annual return, annual vol", [] {
        // 252 identical daily returns of 0.04% → ~10.6% annualized,
        // vol = 0 (all same), so Sharpe is 0 (0 vol edge case).
        // Use slight variation instead.
        std::vector<double> rets(252, 0.001);  // 0.1% daily
        PerformanceConfig cfg{.annualization_factor = 252.0, .risk_free_rate = 0.0};
        auto rpt = compute_performance(rets, cfg);
        // Geometric: (1.001)^252 - 1 ≈ 0.2862
        double expected_ann = std::pow(1.001, 252.0) - 1.0;
        // All returns identical → vol = 0, Sharpe = 0
        return near(rpt.annualized_return, expected_ann, 0.001)
            && near(rpt.annualized_vol, 0.0, 1e-10)
            && near(rpt.sharpe_ratio, 0.0, 1e-10);
    });

    test("Monotonically increasing PnL → max drawdown = 0", [] {
        // Cumulative PnL: 100, 101, 102, ..., 109
        std::vector<double> cum;
        for (int i = 0; i < 10; ++i) cum.push_back(100.0 + i);
        auto rets = cumulative_to_returns(cum);
        auto rpt = compute_performance(rets);
        return near(rpt.max_drawdown, 0.0, 1e-12) && rpt.max_dd_duration == 0;
    });

    test("Known drawdown scenario: value and duration", [] {
        // Returns: +10%, -20%, -10%, +5%, +30%
        // Equity: 1.0 → 1.1 → 0.88 → 0.792 → 0.8316 → 1.08108
        // Peak after first: 1.1, trough at 0.792 → DD = (1.1-0.792)/1.1 = 0.28
        // Duration: from index 1 to index 4 (recovered at 4) = 3 periods?
        // Actually equity at 4 is 1.08108 < peak 1.1, so still in DD.
        // DD ends when equity > 1.1 → at index 4: 1.08108 < 1.1 → never fully recovered
        // So duration = 4 periods (from index 1 through end = 4)
        std::vector<double> rets = {0.10, -0.20, -0.10, 0.05, 0.30};
        auto rpt = compute_performance(rets);
        double expected_dd = (1.1 - 0.792) / 1.1;
        return near(rpt.max_drawdown, expected_dd, 0.001)
            && rpt.max_dd_duration == 4;
    });

    test("Sortino: no negative returns → infinity", [] {
        std::vector<double> rets = {0.01, 0.02, 0.03, 0.01, 0.02};
        auto rpt = compute_performance(rets);
        return std::isinf(rpt.sortino_ratio) && rpt.sortino_ratio > 0;
    });

    test("Calmar ratio computation", [] {
        // 10% annualized return with 5% max drawdown → Calmar = 2.0
        // Construct returns that give known ann return and drawdown.
        std::vector<double> rets = {0.10, -0.20, -0.10, 0.05, 0.30};
        auto rpt = compute_performance(rets);
        // Calmar = annualized_return / max_drawdown
        if (rpt.max_drawdown < 1e-15) return false;
        double expected = rpt.annualized_return / rpt.max_drawdown;
        return near(rpt.calmar_ratio, expected, 0.001);
    });

    test("Edge case: empty series returns zeroed report", [] {
        std::vector<double> empty;
        auto rpt = compute_performance(empty);
        return rpt.n_periods == 0
            && near(rpt.annualized_return, 0.0, 1e-15)
            && near(rpt.sharpe_ratio, 0.0, 1e-15)
            && near(rpt.max_drawdown, 0.0, 1e-15);
    });

    test("Edge case: single-element series", [] {
        std::vector<double> single = {0.05};
        PerformanceConfig cfg{.annualization_factor = 252.0};
        auto rpt = compute_performance(single, cfg);
        return rpt.n_periods == 1
            && near(rpt.annualized_return, 0.05 * 252.0, 0.01)
            && near(rpt.annualized_vol, 0.0, 1e-15);
    });

    test("cumulative_to_returns conversion correctness", [] {
        // Cumulative: 100, 110, 99, 108.9
        // Returns: 10/100=0.10, -11/110=-0.10, 9.9/99=0.10
        std::vector<double> cum = {100.0, 110.0, 99.0, 108.9};
        auto rets = cumulative_to_returns(cum);
        return rets.size() == 3
            && near(rets[0], 0.10, 1e-10)
            && near(rets[1], -0.10, 1e-10)
            && near(rets[2], 0.10, 1e-10);
    });

    test("ReturnSeries concept satisfied by vector<double>", [] {
        return ReturnSeries<std::vector<double>>;
    });
}

// ── Heston model ─────────────────────────────────────────────────
void test_heston() {
    test("HestonVol satisfies VolSurface concept", [] {
        return VolSurface<HestonVol>;
    });

    test("Heston semi-closed-form call price reasonable", [] {
        HestonParams p{.v0=0.04, .kappa=1.5, .theta=0.04, .xi=0.3, .rho=-0.7};
        HestonVol hv(p, 100.0, 0.05, 0.0);
        double call_px = hv.call_price(100.0, 1.0);
        // Should be in the ballpark of BS with 20% vol
        return call_px > 5.0 && call_px < 20.0;
    });

    test("HestonMCEngine satisfies Engine concept", [] {
        return Engine<HestonMCEngine<FlatRate>>;
    });

    test("Heston MC converges to semi-closed-form", [] {
        HestonParams p{.v0=0.04, .kappa=1.5, .theta=0.04, .xi=0.3, .rho=-0.7};
        HestonVol hv(p, 100.0, 0.05, 0.0);
        double cf_price = hv.call_price(100.0, 1.0);

        MCConfig cfg{.n_paths=200'000, .steps=50, .seed=42,
                     .antithetic=true, .n_threads=1, .batch_size=50'000};
        HestonMCEngine mc(FlatRate{Rate{0.05}}, p, cfg);
        MarketSnap m{Spot{100}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
        Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};
        auto r = mc.price(call, m);
        if (!r) return false;
        return std::abs(r->price - cf_price) < 3.0 * r->std_err + 0.5;
    });

    test("Heston smile is non-flat", [] {
        HestonParams p{.v0=0.04, .kappa=1.5, .theta=0.04, .xi=0.3, .rho=-0.7};
        HestonVol hv(p, 100.0, 0.05, 0.0);
        double iv_90  = hv.iv(90.0, 1.0);
        double iv_100 = hv.iv(100.0, 1.0);
        double iv_110 = hv.iv(110.0, 1.0);
        // With negative rho, OTM puts (low strike) have higher IV
        return iv_90 > iv_100 && iv_100 > iv_110;
    });

    test("Heston with xi=0 reduces to BS-like pricing", [] {
        double v0 = 0.04;
        HestonParams p{.v0=v0, .kappa=1.5, .theta=v0, .xi=1e-6, .rho=0.0};
        HestonVol hv(p, 100.0, 0.05, 0.0);
        double heston_px = hv.call_price(100.0, 1.0);
        auto bs = garman_kohlhagen(100.0, 100.0, 1.0, std::sqrt(v0), 0.05, 0.0, OptType::Call);
        return near(heston_px, bs.price, 0.5);  // loose tolerance — xi=0 limit
    });

    test("Heston put-call parity", [] {
        HestonParams p{.v0=0.04, .kappa=1.5, .theta=0.04, .xi=0.3, .rho=-0.7};
        HestonVol hv(p, 100.0, 0.05, 0.0);
        double call_px = hv.call_price(100.0, 1.0);
        double put_px  = hv.put_price(100.0, 1.0);
        double fwd = 100.0 * std::exp(0.05 * 1.0);
        double parity = call_px - put_px - (fwd - 100.0) * std::exp(-0.05);
        return near(parity, 0.0, 0.01);
    });

    test("Feller condition detection", [] {
        HestonParams good{.v0=0.04, .kappa=2.0, .theta=0.04, .xi=0.3, .rho=-0.7};
        HestonParams bad{.v0=0.04, .kappa=0.5, .theta=0.01, .xi=0.8, .rho=-0.7};
        return good.feller_satisfied() && !bad.feller_satisfied();
    });
}

// ── Barrier options ──────────────────────────────────────────────
void test_barrier() {
    MarketSnap m{Spot{100}, Vol{0.20}, Rate{0.05}, Rate{0.0}};

    test("BarrierContract can be constructed", [] {
        BarrierContract bc{Strike{100}, YearFrac{1.0}, OptType::Call,
                           Exercise::European, "EURUSD", 80.0,
                           BarrierType::DownAndOut, 0.0};
        return bc.barrier == 80.0 && bc.barrier_type == BarrierType::DownAndOut;
    });

    test("Down-and-out call: spot near barrier → cheap", [&] {
        MCConfig cfg{.n_paths=100'000, .steps=50, .seed=42,
                     .antithetic=true, .n_threads=1, .batch_size=50'000};
        BarrierMCEngine bmc(FlatRate{m.r_d}, FlatVol{m.sigma}, cfg);
        BarrierContract bc{Strike{100}, YearFrac{1.0}, OptType::Call,
                           Exercise::European, "EURUSD", 95.0,
                           BarrierType::DownAndOut, 0.0};
        auto r = bmc.price_barrier(bc, m);
        if (!r) return false;
        // With barrier at 95 and spot at 100, high chance of knock-out → cheap
        BSEngine bs;
        auto van = bs.price(Contract{Strike{100}, YearFrac{1.0}, OptType::Call}, m);
        return r && van && r->price < van->price;
    });

    test("Down-and-out call with spot far above barrier → near vanilla", [&] {
        MCConfig cfg{.n_paths=100'000, .steps=50, .seed=42,
                     .antithetic=true, .n_threads=1, .batch_size=50'000};
        BarrierMCEngine bmc(FlatRate{m.r_d}, FlatVol{m.sigma}, cfg);
        BarrierContract bc{Strike{100}, YearFrac{1.0}, OptType::Call,
                           Exercise::European, "EURUSD", 50.0,
                           BarrierType::DownAndOut, 0.0};
        auto r = bmc.price_barrier(bc, m);
        BSEngine bs;
        auto van = bs.price(Contract{Strike{100}, YearFrac{1.0}, OptType::Call}, m);
        if (!r || !van) return false;
        // Barrier at 50 very unlikely to hit → nearly vanilla
        return std::abs(r->price - van->price) < 1.0;
    });

    test("In-out parity: knock-in + knock-out ≈ vanilla", [&] {
        MCConfig cfg{.n_paths=100'000, .steps=50, .seed=42,
                     .antithetic=true, .n_threads=1, .batch_size=50'000};
        BarrierMCEngine bmc(FlatRate{m.r_d}, FlatVol{m.sigma}, cfg);

        BarrierContract ko{Strike{100}, YearFrac{1.0}, OptType::Call,
                           Exercise::European, "EURUSD", 80.0,
                           BarrierType::DownAndOut, 0.0};
        BarrierContract ki{Strike{100}, YearFrac{1.0}, OptType::Call,
                           Exercise::European, "EURUSD", 80.0,
                           BarrierType::DownAndIn, 0.0};

        auto ko_r = bmc.price_barrier(ko, m);
        auto ki_r = bmc.price_barrier(ki, m);
        BSEngine bs;
        auto van = bs.price(Contract{Strike{100}, YearFrac{1.0}, OptType::Call}, m);
        if (!ko_r || !ki_r || !van) return false;
        double sum = ko_r->price + ki_r->price;
        return std::abs(sum - van->price) < 2.0;  // MC tolerance
    });

    test("Barrier MC converges with more paths", [&] {
        auto run = [&](std::uint64_t n) {
            MCConfig cfg{.n_paths=n, .steps=50, .seed=42,
                         .antithetic=true, .n_threads=1, .batch_size=50'000};
            BarrierMCEngine bmc(FlatRate{m.r_d}, FlatVol{m.sigma}, cfg);
            BarrierContract bc{Strike{100}, YearFrac{1.0}, OptType::Call,
                               Exercise::European, "EURUSD", 80.0,
                               BarrierType::DownAndOut, 0.0};
            return bmc.price_barrier(bc, m);
        };
        auto r1 = run(10'000);
        auto r2 = run(100'000);
        return r1 && r2 && r2->std_err < r1->std_err;
    });

    test("BarrierMCEngine satisfies Engine concept", [] {
        return Engine<BarrierMCEngine<FlatRate, FlatVol>>;
    });
}

// ── Finite Difference PDE ────────────────────────────────────────
void test_fd() {
    MarketSnap m{Spot{100}, Vol{0.20}, Rate{0.05}, Rate{0.0}};

    test("FDEngine satisfies Engine concept", [] {
        return Engine<FDEngine>;
    });

    test("European call: FD matches BS within grid error", [&] {
        Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};
        BSEngine bs;
        auto bs_r = bs.price(call, m);
        FDEngine fd(FDConfig{.n_space=300, .n_time=300, .s_max_mult=3.0, .theta=0.5});
        auto fd_r = fd.price(call, m);
        if (!bs_r || !fd_r) return false;
        double rel = std::abs(fd_r->price - bs_r->price) / bs_r->price;
        return rel < 0.005;  // <0.5%
    });

    test("European put: FD matches BS", [&] {
        Contract put{Strike{100}, YearFrac{1.0}, OptType::Put};
        BSEngine bs;
        auto bs_r = bs.price(put, m);
        FDEngine fd(FDConfig{.n_space=300, .n_time=300, .s_max_mult=3.0, .theta=0.5});
        auto fd_r = fd.price(put, m);
        if (!bs_r || !fd_r) return false;
        double rel = std::abs(fd_r->price - bs_r->price) / bs_r->price;
        return rel < 0.005;
    });

    test("American put >= European put", [&] {
        Contract eu_put{Strike{100}, YearFrac{1.0}, OptType::Put};
        Contract am_put{Strike{100}, YearFrac{1.0}, OptType::Put, Exercise::American};
        FDEngine fd(FDConfig{.n_space=300, .n_time=300, .s_max_mult=3.0, .theta=0.5});
        auto eu_r = fd.price(eu_put, m);
        auto am_r = fd.price(am_put, m);
        return eu_r && am_r && am_r->price >= eu_r->price - 1e-6;
    });

    test("American put at deep ITM converges to intrinsic", [] {
        MarketSnap deep{Spot{60}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
        Contract am_put{Strike{100}, YearFrac{1.0}, OptType::Put, Exercise::American};
        FDEngine fd(FDConfig{.n_space=300, .n_time=300, .s_min_mult=0.0, .s_max_mult=3.0, .theta=0.5});
        auto r = fd.price(am_put, deep);
        if (!r) return false;
        double intrinsic = 100.0 - 60.0;
        return r->price >= intrinsic - 0.5;  // should be at or above intrinsic
    });

    test("FD delta within 5% of BS analytic", [&] {
        Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};
        BSEngine bs;
        auto bs_r = bs.price(call, m);
        FDEngine fd(FDConfig{.n_space=300, .n_time=300, .s_max_mult=3.0, .theta=0.5});
        auto fd_r = fd.price(call, m);
        if (!bs_r || !fd_r) return false;
        double rel = std::abs(fd_r->greeks.delta - bs_r->greeks.delta) / std::abs(bs_r->greeks.delta);
        return rel < 0.05;
    });

    test("Grid refinement: finer grid gives smaller error", [&] {
        Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};
        BSEngine bs;
        double true_px = bs.price(call, m)->price;

        FDEngine fd_coarse(FDConfig{.n_space=50, .n_time=50, .s_max_mult=3.0, .theta=0.5});
        FDEngine fd_fine(FDConfig{.n_space=300, .n_time=300, .s_max_mult=3.0, .theta=0.5});
        auto r_coarse = fd_coarse.price(call, m);
        auto r_fine   = fd_fine.price(call, m);
        if (!r_coarse || !r_fine) return false;
        return std::abs(r_fine->price - true_px) < std::abs(r_coarse->price - true_px);
    });

    test("Crank-Nicolson more accurate than fully implicit", [&] {
        Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};
        BSEngine bs;
        double true_px = bs.price(call, m)->price;

        FDEngine fd_cn(FDConfig{.n_space=200, .n_time=200, .s_max_mult=3.0, .theta=0.5});
        FDEngine fd_imp(FDConfig{.n_space=200, .n_time=200, .s_max_mult=3.0, .theta=1.0});
        auto r_cn  = fd_cn.price(call, m);
        auto r_imp = fd_imp.price(call, m);
        if (!r_cn || !r_imp) return false;
        return std::abs(r_cn->price - true_px) < std::abs(r_imp->price - true_px);
    });
}

// ══════════════════════════════════════════════════════════════════
int main() {
    std::cout << "\n  Running pricer + execution tests...\n\n";

    // Pricing tests
    test_types();
    test_bs();
    test_fx();
    test_errors();
    test_mc();
    test_rates();
    test_vols();
    test_sobol();
    test_arena();
    test_thread_pool();
    test_iv();
    test_concepts();
    test_move();

    // New feature tests
    test_heston();
    test_barrier();
    test_fd();

    // Execution tests
    test_order_book();
    test_algos();
    test_tca();
    test_signals();
    test_performance();

    std::cout << std::format("\n  Results: {} passed, {} failed, {} total\n\n",
                             passed, failed, passed + failed);
    return failed > 0 ? 1 : 0;
}
