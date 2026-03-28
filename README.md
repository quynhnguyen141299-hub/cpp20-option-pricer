# C++20 FX Option Pricer & Execution Engine

A production-quality library that prices FX options and simulates algorithmic order execution with transaction cost analysis. Built in modern C++ (C++23 / GCC 14+), header-only, zero warnings, 92 tests.

---

## What it does

### Pricing

Values European FX options with multiple engines and cross-validates them:

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
- **Performance analytics** — computes annualized return, volatility, Sharpe ratio, Sortino ratio, Calmar ratio, max drawdown (value and duration) from any PnL or returns series. Concept-constrained, works with any range of doubles.

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

Negative ρ → downside skew (OTM put IV > ATM > OTM call)
```

### Barrier options — Down-and-Out Call

```
 Barrier   DO Price         SE    % Vanilla
  --------------------------------------------
    60.0    10.4679     0.0330       100.2%
    80.0    10.3303     0.0329        98.8%
    90.0     8.6510     0.0326        82.8%
    95.0     5.6430     0.0292        54.0%
```

### FD PDE solver — European & American put

```
            Method      Price      Delta      Gamma
  --------------------------------------------------
       BS European     5.5735    -0.3632   0.018762
  FD European (CN)     5.5727    -0.3618   0.018725
  FD American (CN)     6.0874    -0.4092   0.022923

  FD vs BS error:         0.000831 (0.015%)
  Early exercise premium: 0.514657
```

### Algo comparison (10M EURUSD buy, 1 hour)

```
Algorithm           Fills  Avg Price  Arriv bps   VWAP bps     Impact
-----------------------------------------------------------------------
TWAP                   19    1.08499      -0.10      +0.55      +0.55
VWAP                   20    1.08499      -0.10      +0.46      +0.46
IS-Almgren             19    1.08504      +0.38      +1.01      +1.01
```

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
| **Thomas algorithm** | O(N) tridiagonal solve — the 1D Kronecker specialisation of the FD PDE operator. |

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

src/main.cpp               18 demos (9 pricing + 4 new features + 5 execution)
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


