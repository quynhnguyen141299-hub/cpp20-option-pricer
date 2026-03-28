/// @file main.cpp
/// Demonstration: BS vs MC with convergence, variance reduction, threading.

#include "pricer/pricer.hpp"
#include <iostream>
#include <format>
#include <cmath>
#include <chrono>
#include <vector>
#include <random>

using namespace pricer;

void header(const char* title) {
    std::cout << "\n" << std::string(74, '=')
              << "\n  " << title
              << "\n" << std::string(74, '=') << "\n";
}

// ──────────────────────────────────────────────────────────────────
// Demo 1: FX option — BS analytic vs MC, put-call parity
// ──────────────────────────────────────────────────────────────────
void demo_fx() {
    header("EURUSD 3M 1.0900 — Garman-Kohlhagen vs Monte Carlo");

    MarketSnap mkt{Spot{1.0850}, Vol{0.0750}, Rate{0.0435}, Rate{0.0250}};
    Contract call{Strike{1.0900}, YearFrac{0.25}, OptType::Call, Exercise::European, "EURUSD"};
    Contract put {Strike{1.0900}, YearFrac{0.25}, OptType::Put,  Exercise::European, "EURUSD"};

    std::cout << std::format("  S={:.4f} K={:.4f} T={:.2f} σ={:.2f}% r_d={:.2f}% r_f={:.2f}%\n",
        mkt.S.v, call.K.v, call.T.v, mkt.sigma.v*100, mkt.r_d.v*100, mkt.r_f.v*100);

    BSEngine bs;
    auto bc = bs.price(call, mkt);
    auto bp = bs.price(put, mkt);
    if (bc) std::cout << "  BS Call: " << bc->summary() << "\n";
    if (bp) std::cout << "  BS Put:  " << bp->summary() << "\n";

    // Put-call parity
    if (bc && bp) {
        double lhs = bc->price - bp->price;
        double rhs = mkt.S.v * std::exp(-mkt.r_f.v * call.T.v)
                   - call.K.v * std::exp(-mkt.r_d.v * call.T.v);
        std::cout << std::format("  Parity: C-P={:.8f}  S·df_f-K·df_d={:.8f}  err={:.2e}\n",
                                  lhs, rhs, std::abs(lhs - rhs));
    }

    // MC with control variate + antithetic
    MCConfig cfg{.n_paths = 500'000, .seed = 42, .antithetic = true,
                 .control_variate = true, .n_threads = 1, .batch_size = 50'000};
    MCEngine mc(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg);
    auto mc_c = mc.price(call, mkt);
    auto mc_p = mc.price(put, mkt);
    if (mc_c) std::cout << "  MC Call: " << mc_c->summary() << "\n";
    if (mc_p) std::cout << "  MC Put:  " << mc_p->summary() << "\n";

    if (bc && mc_c) {
        double diff = std::abs(bc->price - mc_c->price);
        std::cout << std::format("  BS-MC diff: {:.2e} ({:.1f} SE)\n",
            diff, mc_c->std_err > 0 ? diff / mc_c->std_err : 0.0);
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 2: Variance reduction comparison
// ──────────────────────────────────────────────────────────────────
void demo_variance_reduction() {
    header("Variance Reduction: Raw vs Antithetic vs Control Variate vs Both");

    MarketSnap mkt{Spot{100.0}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    Contract call{Strike{100.0}, YearFrac{1.0}, OptType::Call};

    BSEngine bs;
    auto analytic = bs.price(call, mkt);
    double true_px = analytic ? analytic->price : 0;
    std::cout << std::format("  Analytic: {:.6f}\n\n", true_px);

    struct VRTest { const char* label; bool at; bool cv; };
    VRTest tests[] = {
        {"Raw (Sobol only)",           false, false},
        {"+ Antithetic",               true,  false},
        {"+ Control Variate",          false, true},
        {"+ Antithetic + CV (full)",   true,  true},
    };

    std::cout << std::format("  {:<30s} {:>10s} {:>10s} {:>10s} {:>10s}\n",
                             "Method", "Price", "Std Err", "|Error|", "µs");
    std::cout << "  " << std::string(72, '-') << "\n";

    for (auto& [label, at, cv] : tests) {
        MCConfig cfg{.n_paths = 200'000, .seed = 42, .antithetic = at,
                     .control_variate = cv, .n_threads = 1, .batch_size = 50'000};
        MCEngine mc(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg);
        auto r = mc.price(call, mkt);
        if (r) {
            std::cout << std::format("  {:<30s} {:>10.6f} {:>10.6f} {:>10.6f} {:>10.0f}\n",
                label, r->price, r->std_err, std::abs(r->price - true_px), r->elapsed_us);
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 3: Convergence — Sobol vs Pseudo-random
// ──────────────────────────────────────────────────────────────────
void demo_convergence() {
    header("Convergence: Sobol QMC vs O(1/√N) pseudo-random baseline");

    MarketSnap mkt{Spot{100.0}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    Contract call{Strike{100.0}, YearFrac{1.0}, OptType::Call};
    BSEngine bs;
    double true_px = bs.price(call, mkt)->price;

    std::cout << std::format("  Analytic: {:.6f}\n\n", true_px);
    std::cout << std::format("  {:>10s} {:>12s} {:>12s} {:>12s} {:>10s}\n",
                             "Paths", "MC Price", "Std Err", "|Error|", "µs");
    std::cout << "  " << std::string(60, '-') << "\n";

    for (auto n : {1'000UL, 5'000UL, 10'000UL, 50'000UL, 100'000UL, 500'000UL, 1'000'000UL}) {
        MCConfig cfg{.n_paths = n, .seed = 42, .antithetic = true,
                     .control_variate = true, .n_threads = 1, .batch_size = 50'000};
        MCEngine mc(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg);
        auto r = mc.price(call, mkt);
        if (r) {
            std::cout << std::format("  {:>10d} {:>12.6f} {:>12.6f} {:>12.6f} {:>10.0f}\n",
                n, r->price, r->std_err, std::abs(r->price - true_px), r->elapsed_us);
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 4: Pathwise Greeks vs Analytic
// ──────────────────────────────────────────────────────────────────
void demo_greeks() {
    header("Pathwise Greeks vs Closed-Form");

    MarketSnap mkt{Spot{100.0}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    Contract call{Strike{100.0}, YearFrac{1.0}, OptType::Call};

    BSEngine bs;
    auto an = bs.price(call, mkt);

    MCConfig cfg{.n_paths = 1'000'000, .seed = 42, .antithetic = true,
                 .control_variate = true, .n_threads = 1, .batch_size = 100'000};
    MCEngine mc(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg);
    auto mcr = mc.price(call, mkt);

    if (an && mcr) {
        std::cout << std::format("  {:>12s} {:>12s} {:>12s} {:>12s}\n",
                                 "Greek", "Analytic", "MC Pathwise", "Error");
        std::cout << "  " << std::string(52, '-') << "\n";
        std::cout << std::format("  {:>12s} {:>12.6f} {:>12.6f} {:>12.6f}\n",
            "Delta", an->greeks.delta, mcr->greeks.delta,
            std::abs(an->greeks.delta - mcr->greeks.delta));
        std::cout << std::format("  {:>12s} {:>12.6f} {:>12.6f} {:>12.6f}\n",
            "Gamma", an->greeks.gamma, mcr->greeks.gamma,
            std::abs(an->greeks.gamma - mcr->greeks.gamma));
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 5: Threading speedup
// ──────────────────────────────────────────────────────────────────
void demo_threading() {
    header("Threading Speedup (1M paths)");

    MarketSnap mkt{Spot{100.0}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    Contract call{Strike{100.0}, YearFrac{1.0}, OptType::Call};

    std::cout << std::format("  {:>8s} {:>12s} {:>12s} {:>12s}\n",
                             "Threads", "Price", "µs", "Speedup");
    std::cout << "  " << std::string(48, '-') << "\n";

    double base_us = 0;
    for (unsigned nt : {1u, 2u, 4u}) {
        MCConfig cfg{.n_paths = 1'000'000, .seed = 42, .antithetic = true,
                     .control_variate = false, .n_threads = nt, .batch_size = 100'000};
        MCEngine mc(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg);
        auto r = mc.price(call, mkt);
        if (r) {
            if (nt == 1) base_us = r->elapsed_us;
            std::cout << std::format("  {:>8d} {:>12.6f} {:>12.0f} {:>12.2f}x\n",
                nt, r->price, r->elapsed_us,
                base_us / std::max(r->elapsed_us, 1.0));
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 6: Early termination via coroutine
// ──────────────────────────────────────────────────────────────────
void demo_early_stop() {
    header("Coroutine Early Termination (target SE = 0.05)");

    MarketSnap mkt{Spot{100.0}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    Contract call{Strike{100.0}, YearFrac{1.0}, OptType::Call};

    // Use NO control variate so SE is meaningful, and a generous target
    // so the coroutine terminates well before exhausting the 5M budget.
    MCConfig cfg{.n_paths = 5'000'000, .seed = 42, .antithetic = true,
                 .control_variate = false, .target_se = 0.05,
                 .n_threads = 1, .batch_size = 5'000};
    MCEngine mc(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg);
    auto r = mc.price(call, mkt);

    if (r) {
        std::cout << "  With early stop:    " << r->summary() << "\n";
    }

    // Compare: without early stop, same budget
    MCConfig cfg2 = cfg;
    cfg2.target_se = 0;
    MCEngine mc2(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg2);
    auto r2 = mc2.price(call, mkt);
    if (r2) {
        std::cout << "  Without early stop: " << r2->summary() << "\n";
    }

    if (r && r2) {
        std::cout << std::format("  Speedup: {:.1f}x (stopped early once SE < 0.05)\n",
            r2->elapsed_us / std::max(r->elapsed_us, 1.0));
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 7: Error handling — std::expected monadic chaining
// ──────────────────────────────────────────────────────────────────
void demo_errors() {
    header("Error Handling: std::expected");

    BSEngine bs;
    Contract call{Strike{100.0}, YearFrac{1.0}, OptType::Call};

    // Negative spot
    MarketSnap bad{Spot{-100}, Vol{0.2}, Rate{0.05}, Rate{0.0}};
    auto r = bs.price(call, bad);
    if (!r) std::cout << "  Neg spot:   " << r.error().str() << "\n";

    // Zero vol
    MarketSnap bad2{Spot{100}, Vol{0.0}, Rate{0.05}, Rate{0.0}};
    auto r2 = bs.price(call, bad2);
    if (!r2) std::cout << "  Zero vol:   " << r2.error().str() << "\n";

    // Monadic chaining
    MarketSnap good{Spot{100}, Vol{0.2}, Rate{0.05}, Rate{0.0}};
    auto chain = bs.price(call, good)
        .transform([](PricingResult pr) { return pr.price; })
        .and_then([](double p) -> Result<std::string> {
            if (p < 0) return std::unexpected(PricingError::at(Err::Overflow, "neg price"));
            return std::format("Validated: {:.4f}", p);
        });
    if (chain) std::cout << "  Chain OK:   " << *chain << "\n";
}

// ──────────────────────────────────────────────────────────────────
// Demo 8: Implied vol solver
// ──────────────────────────────────────────────────────────────────
void demo_iv() {
    header("Implied Vol Solver (Newton-Raphson)");

    MarketSnap mkt{Spot{1.0850}, Vol{0.0750}, Rate{0.0435}, Rate{0.0250}};
    Contract call{Strike{1.0900}, YearFrac{0.25}, OptType::Call, Exercise::European, "EURUSD"};

    BSEngine bs;
    auto px = bs.price(call, mkt);
    if (px) {
        std::cout << std::format("  Target (from σ=7.50%): {:.8f}\n", px->price);
        MarketSnap guess = mkt;
        guess.sigma = Vol{0.30};
        auto iv = bs.implied_vol(call, guess, px->price);
        if (iv) {
            std::cout << std::format("  Recovered: {:.6f}% (err={:.2e})\n",
                                      *iv * 100, std::abs(*iv - 0.075));
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 9: Arena allocator
// ──────────────────────────────────────────────────────────────────
void demo_arena() {
    header("Arena Allocator");

    Arena arena(1024 * 1024);  // 1 MB
    std::cout << std::format("  Capacity: {} bytes\n", arena.capacity());

    auto s1 = arena.alloc_doubles(1000);
    std::cout << std::format("  Alloc 1000 doubles: {} bytes used, span size={}\n",
                              arena.used(), s1.size());

    auto s2 = arena.alloc_doubles(5000);
    std::cout << std::format("  Alloc 5000 doubles: {} bytes used, span size={}\n",
                              arena.used(), s2.size());

    // Write to span (bounds-safe)
    for (std::size_t i = 0; i < s1.size(); ++i) s1[i] = static_cast<double>(i);
    std::cout << std::format("  s1[999] = {:.0f}\n", s1[999]);

    arena.reset();
    std::cout << std::format("  After reset: {} bytes used\n", arena.used());
}

// ──────────────────────────────────────────────────────────────────
// Demo 10: Heston stochastic vol smile generation
// ──────────────────────────────────────────────────────────────────
void demo_heston() {
    header("Heston Stochastic Vol — Smile Generation");

    HestonParams p{.v0=0.04, .kappa=1.5, .theta=0.04, .xi=0.3, .rho=-0.7};
    double S = 100.0, r_d = 0.05, r_f = 0.0, T = 0.5;
    HestonVol hv(p, S, r_d, r_f);

    std::cout << std::format("  v0={} κ={} θ={} ξ={} ρ={}\n",
        p.v0, p.kappa, p.theta, p.xi, p.rho);
    std::cout << std::format("  Feller condition: {}\n\n",
        p.feller_satisfied() ? "satisfied" : "VIOLATED");

    std::cout << std::format("  {:>8s} {:>10s} {:>12s} {:>12s}\n",
                             "Strike", "Moneyness", "Call Price", "Impl Vol");
    std::cout << "  " << std::string(46, '-') << "\n";

    for (double m : {0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20}) {
        double K = S * m;
        double px = hv.call_price(K, T);
        double iv = hv.iv(K, T);
        std::cout << std::format("  {:>8.2f} {:>10.2f} {:>12.4f} {:>11.2f}%\n",
                                  K, m, px, iv * 100);
    }
    std::cout << "\n  Negative ρ → downside skew (OTM put IV > ATM > OTM call)\n";
}

// ──────────────────────────────────────────────────────────────────
// Demo 11: Heston MC vs semi-closed-form
// ──────────────────────────────────────────────────────────────────
void demo_heston_mc() {
    header("Heston MC (QE Scheme) vs Semi-Closed-Form");

    HestonParams p{.v0=0.04, .kappa=1.5, .theta=0.04, .xi=0.3, .rho=-0.7};
    double S = 100.0, r_d = 0.05, r_f = 0.0;
    HestonVol hv(p, S, r_d, r_f);

    MarketSnap mkt{Spot{S}, Vol{0.20}, Rate{r_d}, Rate{r_f}};
    Contract call{Strike{100}, YearFrac{1.0}, OptType::Call};

    double cf_price = hv.call_price(100.0, 1.0);
    std::cout << std::format("  Semi-closed-form: {:.6f}\n\n", cf_price);

    std::cout << std::format("  {:>10s} {:>8s} {:>12s} {:>10s} {:>10s}\n",
                             "Paths", "Steps", "MC Price", "SE", "|Error|");
    std::cout << "  " << std::string(54, '-') << "\n";

    for (auto [n, s] : std::vector<std::pair<uint64_t,uint32_t>>{
        {50'000, 50}, {100'000, 50}, {200'000, 50}, {500'000, 50}
    }) {
        MCConfig cfg{.n_paths=n, .steps=s, .seed=42, .n_threads=1};
        HestonMCEngine mc(FlatRate{Rate{r_d}}, p, cfg);
        auto r = mc.price(call, mkt);
        if (r) {
            std::cout << std::format("  {:>10d} {:>8d} {:>12.6f} {:>10.4f} {:>10.4f}\n",
                n, s, r->price, r->std_err, std::abs(r->price - cf_price));
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 12: Barrier option pricing
// ──────────────────────────────────────────────────────────────────
void demo_barrier() {
    header("Barrier Options — Down-and-Out Call with Brownian Bridge");

    MarketSnap mkt{Spot{100.0}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    BSEngine bs;
    auto vanilla = bs.price(Contract{Strike{100}, YearFrac{1.0}, OptType::Call}, mkt);
    if (vanilla) {
        std::cout << std::format("  Vanilla call: {:.4f}\n\n", vanilla->price);
    }

    MCConfig cfg{.n_paths=200'000, .steps=100, .seed=42,
                 .antithetic=false, .n_threads=1};
    BarrierMCEngine bmc(FlatRate{mkt.r_d}, FlatVol{mkt.sigma}, cfg);

    std::cout << std::format("  {:>8s} {:>10s} {:>10s} {:>12s}\n",
                             "Barrier", "DO Price", "SE", "% Vanilla");
    std::cout << "  " << std::string(44, '-') << "\n";

    for (double B : {60.0, 70.0, 80.0, 85.0, 90.0, 95.0}) {
        BarrierContract bc{Strike{100}, YearFrac{1.0}, OptType::Call,
                           Exercise::European, "EURUSD", B,
                           BarrierType::DownAndOut, 0.0};
        auto r = bmc.price_barrier(bc, mkt);
        if (r && vanilla) {
            std::cout << std::format("  {:>8.1f} {:>10.4f} {:>10.4f} {:>11.1f}%\n",
                B, r->price, r->std_err, 100.0 * r->price / vanilla->price);
        }
    }

    std::cout << "\n  Barrier closer to spot → more knock-out probability → cheaper\n";
}

// ──────────────────────────────────────────────────────────────────
// Demo 13: Finite Difference PDE solver
// ──────────────────────────────────────────────────────────────────
void demo_fd() {
    header("FD PDE Solver — European & American Put");

    MarketSnap mkt{Spot{100.0}, Vol{0.20}, Rate{0.05}, Rate{0.0}};
    BSEngine bs;

    Contract eu_put{Strike{100}, YearFrac{1.0}, OptType::Put};
    Contract am_put{Strike{100}, YearFrac{1.0}, OptType::Put, Exercise::American};

    auto bs_r = bs.price(eu_put, mkt);
    FDEngine fd(FDConfig{.n_space=300, .n_time=300, .s_max_mult=3.0, .theta=0.5});
    auto fd_eu = fd.price(eu_put, mkt);
    auto fd_am = fd.price(am_put, mkt);

    if (bs_r && fd_eu && fd_am) {
        std::cout << std::format("  {:>20s} {:>10s} {:>10s} {:>10s}\n",
                                 "Method", "Price", "Delta", "Gamma");
        std::cout << "  " << std::string(44, '-') << "\n";
        std::cout << std::format("  {:>20s} {:>10.4f} {:>+10.4f} {:>10.6f}\n",
            "BS European", bs_r->price, bs_r->greeks.delta, bs_r->greeks.gamma);
        std::cout << std::format("  {:>20s} {:>10.4f} {:>+10.4f} {:>10.6f}\n",
            "FD European (CN)", fd_eu->price, fd_eu->greeks.delta, fd_eu->greeks.gamma);
        std::cout << std::format("  {:>20s} {:>10.4f} {:>+10.4f} {:>10.6f}\n",
            "FD American (CN)", fd_am->price, fd_am->greeks.delta, fd_am->greeks.gamma);

        double err = std::abs(fd_eu->price - bs_r->price);
        double premium = fd_am->price - fd_eu->price;
        std::cout << std::format("\n  FD vs BS error:         {:.6f} ({:.3f}%)\n",
            err, 100.0 * err / bs_r->price);
        std::cout << std::format("  Early exercise premium: {:.6f}\n", premium);
    }

    // Grid refinement study
    std::cout << "\n  Grid refinement study (European call):\n";
    Contract eu_call{Strike{100}, YearFrac{1.0}, OptType::Call};
    double true_px = bs.price(eu_call, mkt)->price;
    std::cout << std::format("  {:>6s} {:>10s} {:>12s} {:>10s}\n",
                             "N×M", "FD Price", "|Error|", "µs");
    std::cout << "  " << std::string(42, '-') << "\n";
    for (int n : {50, 100, 200, 400}) {
        FDEngine fd_n(FDConfig{.n_space=static_cast<uint32_t>(n),
                               .n_time=static_cast<uint32_t>(n),
                               .s_max_mult=3.0, .theta=0.5});
        auto r = fd_n.price(eu_call, mkt);
        if (r) {
            std::cout << std::format("  {:>3d}×{:<3d} {:>10.6f} {:>12.6f} {:>10.0f}\n",
                n, n, r->price, std::abs(r->price - true_px), r->elapsed_us);
        }
    }
}

// ══════════════════════════════════════════════════════════════════
//  EXECUTION DEMOS
// ══════════════════════════════════════════════════════════════════

using namespace pricer::execution;

// ──────────────────────────────────────────────────────────────────
// Demo 10: TWAP vs VWAP vs IS — algo comparison with TCA
// ──────────────────────────────────────────────────────────────────
void demo_algo_comparison() {
    header("Algo Execution: TWAP vs VWAP vs IS (10M EURUSD buy)");

    ExecutionOrder order{
        .total_qty     = 10e6,
        .side          = Side::Ask,
        .start_time_us = 0,
        .end_time_us   = 3600e6,   // 1 hour
        .label         = "EURUSD"
    };

    ImpactModel impact{.eta = 0.2, .daily_vol = 0.007, .daily_volume = 50e9};
    BacktestConfig cfg{.n_ticks = 5'000, .tick_interval_us = 720'000};

    // Run all three algos on identical market conditions (same seed)
    std::cout << std::format("  Order: BUY {:.0f} {} over {:.0f} min\n",
        order.total_qty, order.label,
        (order.end_time_us - order.start_time_us) / 60e6);
    std::cout << std::format("  Impact model: η={:.2f} σ_daily={:.3f}%\n\n",
        impact.eta, impact.daily_vol * 100);

    // Table header
    std::cout << std::format("  {:<14s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s}\n",
        "Algorithm", "Fills", "Avg Price", "Arriv bps", "VWAP bps", "Impact");
    std::cout << "  " << std::string(66, '-') << "\n";

    auto run_algo = [&](auto algo, std::uint64_t seed) {
        OrderBookSim book({}, seed);
        auto result = run_backtest(algo, order, book, impact, cfg);
        auto& t = result.tca;
        std::cout << std::format("  {:<14s} {:>10d} {:>10.5f} {:>+10.2f} {:>+10.2f} {:>+10.2f}\n",
            t.algo_name, t.n_fills, t.avg_fill_price,
            t.arrival_cost_bps, t.vwap_slip_bps, t.impact_cost_bps);
        return result;
    };

    run_algo(TWAPAlgo(order, 20), 42);
    run_algo(VWAPAlgo(order, 20), 42);
    run_algo(ISAlgo(order, 20, 1.0), 42);
}

// ──────────────────────────────────────────────────────────────────
// Demo 11: IS urgency parameter sweep
// ──────────────────────────────────────────────────────────────────
void demo_urgency_sweep() {
    header("IS Urgency Sweep: impact vs drift trade-off");

    ExecutionOrder order{
        .total_qty = 10e6, .side = Side::Ask,
        .start_time_us = 0, .end_time_us = 3600e6, .label = "EURUSD"
    };
    ImpactModel impact{.eta = 0.2, .daily_vol = 0.007, .daily_volume = 50e9};
    BacktestConfig cfg{.n_ticks = 5'000, .tick_interval_us = 720'000};

    std::cout << std::format("  {:<12s} {:>12s} {:>12s} {:>12s}\n",
        "Urgency", "Arriv bps", "Impact bps", "Timing bps");
    std::cout << "  " << std::string(50, '-') << "\n";

    for (double urgency : {0.1, 0.5, 1.0, 2.0, 5.0}) {
        ISAlgo algo(order, 20, urgency);
        OrderBookSim book({}, 42);
        auto result = run_backtest(algo, order, book, impact, cfg);
        auto& t = result.tca;
        std::cout << std::format("  λ={:<10.1f} {:>+12.2f} {:>+12.2f} {:>+12.2f}\n",
            urgency, t.arrival_cost_bps, t.impact_cost_bps, t.timing_risk_bps);
    }

    std::cout << "\n  Low λ → front-loads (high impact, low drift risk)\n"
                 "  High λ → spreads evenly (low impact, high drift risk)\n";
}

// ──────────────────────────────────────────────────────────────────
// Demo 12: Detailed TCA report
// ──────────────────────────────────────────────────────────────────
void demo_tca_report() {
    header("Full TCA Report — VWAP Algo");

    ExecutionOrder order{
        .total_qty = 25e6, .side = Side::Ask,
        .start_time_us = 0, .end_time_us = 3600e6, .label = "EURUSD"
    };
    ImpactModel impact{.eta = 0.3, .daily_vol = 0.007, .daily_volume = 50e9};
    BacktestConfig cfg{.n_ticks = 5'000, .tick_interval_us = 720'000};

    VWAPAlgo algo(order, 20);
    OrderBookSim book({}, 123);
    auto result = run_backtest(algo, order, book, impact, cfg);

    std::cout << "  " << result.tca.summary() << "\n";
}

// ──────────────────────────────────────────────────────────────────
// Demo 13: Alpha signals
// ──────────────────────────────────────────────────────────────────
void demo_signals() {
    header("Execution Signals: Spread + Momentum");

    OrderBookSim book({}, 42);
    SpreadSignal spread_sig(50);
    MomentumSignal mom_sig(30);

    std::cout << std::format("  {:>8s} {:>10s} {:>10s} {:>12s} {:>12s}\n",
        "Tick", "Mid", "Spread", "SpreadSig", "MomentumSig");
    std::cout << "  " << std::string(54, '-') << "\n";

    for (int i = 0; i < 200; ++i) {
        auto tick = book.step(i * 100'000.0);
        double ss = spread_sig.score(tick);
        double ms = mom_sig.score(tick);

        // Print every 20th tick
        if (i % 20 == 0) {
            std::cout << std::format("  {:>8d} {:>10.5f} {:>10.2f} {:>+12.4f} {:>+12.4f}\n",
                i, tick.mid, tick.spread_bps(), ss, ms);
        }
    }
}

// ──────────────────────────────────────────────────────────────────
// Demo 14: Performance analytics on backtest results
// ──────────────────────────────────────────────────────────────────
void demo_performance() {
    header("Performance Analytics");

    // 1. Synthetic daily PnL: positive-drift strategy (mean ~8bp/day, vol ~1.2%)
    std::mt19937 rng(7777);
    std::normal_distribution<double> dist(0.0008, 0.012);

    std::vector<double> daily_returns(252);  // 1 year of daily returns
    for (auto& r : daily_returns) r = dist(rng);

    auto rpt = compute_performance(daily_returns);
    std::cout << "  --- Synthetic strategy (252 days, μ=8bp/day) ---\n";
    std::cout << "  " << format_report(rpt) << "\n\n";

    // 2. Performance on VWAP backtest results — build PnL from fills
    ExecutionOrder order{
        .total_qty = 10e6, .side = Side::Ask,
        .start_time_us = 0, .end_time_us = 3600e6, .label = "EURUSD"
    };
    ImpactModel impact{.eta = 0.2, .daily_vol = 0.007, .daily_volume = 50e9};
    BacktestConfig bt_cfg{.n_ticks = 5'000, .tick_interval_us = 720'000};

    VWAPAlgo algo(order, 20);
    OrderBookSim book({}, 42);
    auto result = run_backtest(algo, order, book, impact, bt_cfg);

    // Convert tick mid-prices into a cumulative equity-like series
    std::vector<double> equity;
    equity.reserve(result.ticks.size());
    double base_mid = result.ticks.front().mid;
    for (const auto& t : result.ticks) {
        equity.push_back(t.mid / base_mid);
    }
    auto tick_returns = cumulative_to_returns(equity);

    // Treat ticks as ~1-second bars → annualize with trading seconds/year
    PerformanceConfig tick_cfg{.annualization_factor = 252.0 * 6.5 * 3600.0};
    auto bt_rpt = compute_performance(tick_returns, tick_cfg);
    std::cout << "  --- VWAP backtest tick-level metrics ---\n";
    std::cout << "  " << format_report(bt_rpt) << "\n";
}

// ──────────────────────────────────────────────────────────────────
// Demo: Machine Epsilon
// ──────────────────────────────────────────────────────────────────
void demo_machine_epsilon() {
    header("Machine Epsilon");

    constexpr double eps = std::numeric_limits<double>::epsilon();   // 2^{-52}
    constexpr float  eps_f = std::numeric_limits<float>::epsilon();  // 2^{-23}

    std::cout << std::format("  double epsilon (2^-52):  {:.6e}\n", eps);
    std::cout << std::format("  float  epsilon (2^-23):  {:.6e}\n", static_cast<double>(eps_f));

    // Demonstrate: 1.0 + eps != 1.0, but 1.0 + eps/2 == 1.0
    const double a = 1.0 + eps;
    const double b = 1.0 + eps / 2.0;
    std::cout << std::format("  1.0 + eps     == 1.0 ?   {}\n", a == 1.0 ? "true" : "false");
    std::cout << std::format("  1.0 + eps/2   == 1.0 ?   {}\n", b == 1.0 ? "true" : "false");

    // Show how it relates to our pricer: put-call parity gap
    MarketSnap mkt{Spot{1.0850}, Vol{0.0750}, Rate{0.0435}, Rate{0.0250}};
    Contract call{Strike{1.0900}, YearFrac{0.25}, OptType::Call, Exercise::European, "EURUSD"};
    Contract put {Strike{1.0900}, YearFrac{0.25}, OptType::Put,  Exercise::European, "EURUSD"};
    BSEngine bs;
    auto c = bs.price(call, mkt);
    auto p = bs.price(put,  mkt);
    if (c && p) {
        double S  = 1.0850;
        double K  = 1.0900;
        double rd = 0.0435, rf = 0.0250, T = 0.25;
        double parity_gap = std::abs((c->price - p->price)
                            - (S * std::exp(-rf * T) - K * std::exp(-rd * T)));
        std::cout << std::format("  Put-call parity gap:     {:.2e}  (≈ machine epsilon)\n", parity_gap);
    }
}

// ══════════════════════════════════════════════════════════════════
int main() {
    std::cout << "\n"
        "  ╔════════════════════════════════════════════════════════════════╗\n"
        "  ║  C++20 FX Option Pricer & Execution Engine                   ║\n"
        "  ║  Pricing: GK · MC · Heston · Barrier · FD PDE               ║\n"
        "  ║  Execution: TWAP · VWAP · IS · TCA · Alpha Signals          ║\n"
        "  ╚════════════════════════════════════════════════════════════════╝\n";

    // ── Pricing demos ──
    demo_fx();
    demo_variance_reduction();
    demo_convergence();
    demo_greeks();
    demo_threading();
    demo_early_stop();
    demo_errors();
    demo_iv();
    demo_arena();

    // ── New feature demos ──
    demo_heston();
    demo_heston_mc();
    demo_barrier();
    demo_fd();

    // ── Execution demos ──
    demo_algo_comparison();
    demo_urgency_sweep();
    demo_tca_report();
    demo_signals();
    demo_performance();

    // ── Numerical foundations ──
    demo_machine_epsilon();

    std::cout << "\n" << std::string(74, '=')
              << "\n  All demos complete.\n"
              << std::string(74, '=') << "\n\n";
}
