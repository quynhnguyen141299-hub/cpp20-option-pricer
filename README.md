# C++23 FX Option Pricer & Execution Engine 


A production-quality library that prices FX options and simulates algorithmic order execution with transaction cost analysis. Built in C++ (C++23 / GCC 14+).

QR Usage: Model experimentation
Compare numerical methods: Monte Carlo vs PDE vs closed-form
Study Convergence and Error + Test variance reduction techniques of MC estimators for option pricing

Also, this app deliver fast pricing for FX options and their real-time Greeks using variance reduction techniques (SOBOL QMC + Antithetic Variates + Control Variates + Combined (A+CV)); when you cannot simulate a million MC paths every time, which costs memory and run-time.

---

## What it does

### Pricing

Prices European, American, and barrier FX options with multiple engines and cross-validates them: i.e. checks that the different engines agree with each other or not (for QR options modelling experimentation / quick execution pricing with real-time greeks). 

For example, the MC price should converge to the GK analytical price, and the FD PDE price should match BS to within 0.015%. If they all agree, you have confidence the implementations are correct.

- **Garman-Kohlhagen (analytical)** — exact price and all Greeks in microseconds, plus a Newton-Raphson implied vol solver.
- **Monte Carlo** — Sobol QMC, S_T control variate (~3x variance reduction), antithetic variates, pathwise delta, likelihood-ratio gamma, coroutine early-stop (816x speedup), jthread pool parallelism, and arena allocation.
- **Heston Stochastic Volatility** — semi-closed-form pricing via the Heston (1993) characteristic function with 64-point Gauss-Legendre quadrature. MC engine uses the QE (Quadratic Exponential) discretisation scheme by Andersen (2008) with Cholesky correlation. Generates realistic volatility smiles.
- **Barrier Options** — Monte Carlo pricing of knock-in/knock-out barrier options with Brownian bridge continuity correction (Beaglehole-Dyer-Jacka). In-out parity for knock-in efficiency.
- **Finite Difference PDE Solver** — 1D Black-Scholes PDE with non-uniform sinh grid, θ-scheme time stepping (Crank-Nicolson), Thomas algorithm for tridiagonal solve, American option support via early exercise constraint, Greeks extraction from the grid.

### Execution

Simulates how a large FX order (e.g. buy 10M EURUSD) gets sliced and executed over time:

- **TWAP** — splits the order into equal slices at regular intervals. Simple, predictable, widely used for low-urgency orders.
- **VWAP** — weights slices by an intraday volume profile (U-shaped: heavier at open and close). Tracks the VWAP benchmark.
- **Implementation Shortfall (Almgren-Chriss)** — front-loads execution to minimise arrival-price shortfall. The urgency parameter λ controls the trade-off between market impact and drift risk.
- **Order book simulator** — generates realistic tick data with configurable spread, depth, and a square-root market impact model.
- **Transaction cost analysis (TCA)** — post-trade report: arrival shortfall, VWAP slippage, spread cost, impact cost, timing risk, participation rate.
- **Alpha signals** — spread and momentum signals that score market conditions in [-1, +1], used to adjust execution aggressiveness in real time.
- **Performance analytics** — computes Annualized Return, Volatility, Sharpe, Sortino, Calmar, Max DD (value and duration) from any PnL or doubles returns series (the C++ perf. analytics component take either a PnL series or a std::vector<double> of returns and compute a standard set of risk/return metrics).

---

## Demo output

### Heston volatility smile (ρ = −0.7)

```
  Strike  Moneyness   Call Price     Impl Vol
  ----------------------------------------------
    80.00       0.80      22.4553       24.03%
    90.00       0.90      13.8369       21.89%
   100.00       1.00       6.8257       19.77%
   110.00       1.10       2.3341       17.74%
   120.00       1.20       0.4705       16.11%

ρ = −0.7 means spot and vol are negatively correlated — when the price drops, volatility spikes (the "fear" effect you see in real markets).
- IV trending down as price rises → downside skew (OTM put IV > ATM > OTM call): realistic in FX and equity markets. BS give flat 19.77% across strikes, wrong and Heston model capture it. The Heston model lets volatility itself be random (stochastic vol), unlike Black-Scholes where volatility is fixed.

Strike	meaning:
lower: 80 (deep ITM call)	Moneyness 0.80 means this strike is 20% below spot. IV being 24.03% means high premium vol charge. 
mid: 100 (ATM)	At-the-money. Implied vol is 19.77% —  baseline. Black-Scholes would give a flat 19.77% across all strikes, which is wrong --- Heston model captures this.
above: 120 (deep OTM call)	20% above spot. Implied vol drops to 16.11%. Upside options are "cheaper" in vol terms.
```

### Barrier options — Down-and-Out Call
A down-and-out call is a regular call that dies (knocks out) if spot ever touches the barrier from above during the option's life.

```
 Barrier   DO Price         SE    % Vanilla
  --------------------------------------------
    60.0    10.4679     0.0330       100.2%
    80.0    10.3303     0.0329        98.8%
    90.0     8.6510     0.0326        82.8%
    95.0     5.6430     0.0292        54.0%

Barrier	Interpretation
60.0 → 100.2% of vanilla	Barrier is so far below spot that it's almost never hit. The DO call is worth essentially the same as a plain vanilla call. The 100.2% (slightly above 100%) is MC sampling noise.
80.0 → 98.8%	Still far away, tiny chance of knockout. You lose only 1.2% of the vanilla value.
90.0 → 82.8%	Barrier is getting close to spot. Meaningful probability of knockout, so the option is worth 17% less.
95.0 → 54.0%	Barrier is very close — spot only needs to drop 5% to kill the option. Worth roughly half of vanilla. This is the discount you get for accepting knockout risk.
SE is the Monte Carlo standard error — all around 0.03, meaning the prices are accurate to about ±0.06 (2 SE).
```
### FD PDE solver — European & American put
This compares three methods for pricing the same put option:

```
            Method      Price      Delta      Gamma
  --------------------------------------------------
       BS European     5.5735    -0.3632   0.018762
  FD European (CN)     5.5727    -0.3618   0.018725
  FD American (CN)     6.0874    -0.4092   0.022923

  FD vs BS error:         0.000831 (0.015%)
  Early exercise premium: 0.514657

Method and What it is
- BS European (5.5735) :	Exact analytical Black-Scholes price. The "truth" for European puts.
- FD European CN (5.5727) :	The finite-difference PDE grid solver using Crank-Nicolson. It gets within 0.015% of the exact answer — this validates that the PDE solver is working correctly.
- FD American CN (6.0874) :	Same PDE solver but with early exercise allowed at every time step.

Key numbers:
- Early exercise premium = 0.5147 — the American put is worth 0.51 more than the European put.
- That's the value of being able to exercise early (useful when the option is deep ITM and you'd rather have cash earning interest).
```

--- Delta: −0.41 vs −0.36 — the American put has a steeper delta because it's more likely to be exercised early, making it behave more like the underlying.
--- Gamma: 0.0229 vs 0.0187 — higher gamma on the American put means its delta is more sensitive near the early-exercise boundary.

### Algo comparison (10M EURUSD buy, 1 hour)

```
Algorithm           Fills  Avg Price  Arriv bps   VWAP bps     Impact bps
--------------------------------------------------------------------------
TWAP                   19    1.08499      -0.10      +0.55      +0.55
VWAP                   20    1.08499      -0.10      +0.46      +0.46
IS-Almgren             19    1.08504      +0.38      +1.01      +1.01
```
Note: impact = η × σ × √(quantity / daily_volume)
So for each algorithm:
- TWAP +0.55 bps — your equal-sized slices moved the price 0.55 bps against you over the hour.
- VWAP +0.46 bps — volume-weighted slicing caused less impact as it trades heavier when liquidity is deeper (open/close).
- IS-Almgren +1.01 bps — front-loading the order (to reduce timing risk) causes more impact because you're hitting book harder early on

The trade-off IS-Almgren: accepting higher (market impact) now to avoid risk of the (price drifting) further away if you wait. 
The urgency parameter (λ) controls that balance.

### Performance analytics (synthetic daily returns, 252 days)

```
Performance Report (252 periods)
  Annualized return:      +0.6038  (+60.38%)
  Annualized vol:          0.2062  (20.62%)
  Sharpe ratio:           +2.3951
  Sortino ratio:          +3.9347
  Max drawdown:            0.0852  (8.52%)
  Max DD duration:             47  periods
  Calmar ratio:           +7.0837
```

### Machine epsilon

```
double epsilon (2^-52):  2.220446e-16
float  epsilon (2^-23):  1.192093e-07
1.0 + eps     == 1.0 ?   false
1.0 + eps/2   == 1.0 ?   true
Put-call parity gap:     1.11e-16  (≈ machine epsilon)
```

Machine epsilon is the smallest `double` such that `1.0 + ε ≠ 1.0`. Our GK pricer's put-call parity residual lands exactly at this hardware limit — proof the analytical engine is numerically exact.

### Pricing benchmarks

```
Put-call parity error:           1.11e-16 (machine epsilon)
MC vs BS price difference:       6.13e-07 (< 1 standard error)
Pathwise delta error:            0.000006
Coroutine early-stop:            816× faster (90K paths vs 5M budget)
IV solver round-trip error:      4.16e-17
```

---

| Feature | The problem it solves |
|---|---|
| **Concepts** (`RateModel`, `VolSurface`, `Engine`, `ExecutionAlgo`, `Signal`) | Swap implementations (FlatRate↔PiecewiseCurve, TWAP↔IS) without virtual dispatch. Compile-time constraints. |
| **`std::expected`** | Zero exceptions on the hot path. Monadic error propagation. |
| **Coroutines** (`Generator<BatchResult>`) | MC loop lazily produces batches; consumer stops when precision is sufficient. |
| **`std::jthread` + `stop_token`** | RAII thread lifetime, cooperative cancellation. |
| **Move-only types** | SobolEngine (no sequence correlation), Arena (no double-free), OrderBookSim (no RNG duplication). |
| **Arena allocator** | Zero malloc in the MC inner loop. 64-byte cache-line alignment. |
| **Phantom types** (`Spot`, `Strike`, `Vol`, `Rate`) | Prevents parameter transposition bugs at zero runtime cost. |
| **`std::source_location`** | Every error captures file:line at compile time. |
| **Ranges** | Clean algorithm calls on rate curves and vol term structures. |
| **Designated initializers** | Self-documenting config: `MCConfig{.n_paths=500'000, .antithetic=true}`. |
| **Heston semi-closed-form** | Industry-standard stochastic vol model with characteristic function pricing. |
| **QE discretisation** | Andersen (2008) variance scheme — correct non-negativity, second-order accuracy. |
| **Brownian bridge** | Barrier continuity correction — captures between-step crossings, not just discrete monitoring. |
| **Thomas algorithm** | O(N) tridiagonal solve — the 1D specialisation of the FD PDE operator. |

---

## Project structure

```
include/pricer/
├── core/
│   ├── types.hpp          Phantom-typed doubles, Contract, MarketSnap, Greeks, MCConfig
│   ├── errors.hpp         std::expected<T, PricingError> with source_location
│   ├── concepts.hpp       RateModel, VolSurface, Engine concepts
│   ├── sobol.hpp          Gray-code Sobol QMC + Beasley-Springer-Moro inverse normal
│   ├── arena.hpp          Monotonic arena allocator (RAII, move-only)
│   └── thread_pool.hpp    std::jthread pool with stop_token
├── models/
│   ├── rates.hpp          FlatRate, PiecewiseCurve, FXRatePair<D,F>
│   ├── vol.hpp            FlatVol, TermStructureVol, StickyStrikeVol
│   ├── heston.hpp         HestonParams, HestonVol (semi-closed-form + IV)
│   └── barrier.hpp        BarrierType enum, BarrierContract
├── engines/
│   ├── black_scholes.hpp  GK analytic pricer + Greeks + IV solver
│   ├── mc.hpp             Full MC engine (all variance reduction + threading)
│   ├── heston_mc.hpp      Heston MC with QE discretisation (Andersen 2008)
│   ├── barrier_mc.hpp     Barrier MC with Brownian bridge correction
│   └── fd_pde.hpp         1D FD PDE solver (Thomas, sinh grid, American)
├── execution/
│   ├── order_book.hpp     LOB simulator with market impact model
│   ├── algo.hpp           TWAP, VWAP, IS (Almgren-Chriss) execution algos
│   ├── tca.hpp            Transaction cost analysis engine
│   ├── signal.hpp         Alpha signals (spread, momentum, composite)
│   ├── backtest.hpp       Backtest harness: tick generation → algo → fills → TCA
│   └── performance.hpp    Performance analytics (Sharpe, Sortino, Calmar, drawdown)
└── pricer.hpp             Umbrella header

python/
├── bindings.cpp           pybind11 bindings for all key types and engines
└── demo.py                Python demo script

src/main.cpp               19 demos (9 pricing + 4 new features + 5 execution + 1 numerical)
tests/test_pricer.cpp      92 unit tests
CMakeLists.txt
```

---

## Build

```bash
# GCC 14+ (direct)
g++ -std=c++23 -fcoroutines -O2 -Wall -Wextra -Wpedantic \
    -I include -pthread src/main.cpp -o fx_pricer

# Run tests
g++ -std=c++23 -fcoroutines -O2 -Wall -Wextra -Wpedantic \
    -I include -pthread tests/test_pricer.cpp -o test_pricer
./test_pricer

# CMake
cmake -B build && cmake --build build
ctest --test-dir build

# Python bindings (optional, requires pybind11)
cmake -DBUILD_PYTHON=ON -B build && cmake --build build
```

---


