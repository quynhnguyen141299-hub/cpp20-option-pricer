# C++23 FX Option Pricer & Execution Engine 


A reusable C++ library that helps FX traders prices FX options and FX QRs simulates Algo. strategies execution perf. of large FX orders using TWAP, VWAP, and Almgren-Chriss strategies, with order book simulation, square-root market impact, TCA and real-time alpha signals (spread and momentum). 

Built in C++23 (GCC 14+), header-only, zero warnings, 92 tests.

QR Usage: Model experimentation
- Compare numerical methods: Monte Carlo vs PDE vs closed-form
- Study Convergence and Error + Test variance reduction techniques of MC estimators for option pricing

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
  SpreadSignal
  Keeps a std::deque of the last 100 spread observations (configurable via lookback)
  On each tick: score = (avg - current_spread) / avg, clamped to [-1, +1]
  So it's a simple arithmetic mean of the last 100 ticks — not EWMA, not time-bucketed
  
  MomentumSignal
  Keeps a std::deque of the last 50 mid-prices
  Computes return = (latest_mid - oldest_mid) / oldest_mid
  Normalizes by dividing by 0.01 (assumes daily returns rarely exceed 1%), clamped to [-1, +1]
  
  CompositeSignal
  Weighted average of any mix of signals via std::function type erasure
  
  Normalizes by total absolute weight, clamped to [-1, +1]

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

## Getting started

### Prerequisites

You need **GCC 14+** (for C++23 support). No other dependencies — the library is header-only.

### Windows (via WSL)

**1. Install WSL and Ubuntu**

Open PowerShell as Administrator:

```powershell
wsl --install
```

Restart your PC. On reboot, Ubuntu will open — create a username and password.

**2. Install GCC 14**

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y gcc-14 g++-14 git
```

**3. Clone and build**

```bash
git clone https://github.com/quynhnguyen141299-hub/cpp20-option-pricer.git
cd cpp20-option-pricer

g++-14 -std=c++23 -fcoroutines -O2 -Wall -Wextra -Wpedantic \
    -I include -pthread src/main.cpp -o fx_pricer

./fx_pricer
```

### Linux (native)

Same as steps 2–3 above. Skip the WSL part.

### Mac

```bash
brew install gcc@14
git clone https://github.com/quynhnguyen141299-hub/cpp20-option-pricer.git
cd cpp20-option-pricer

g++-14 -std=c++23 -fcoroutines -O2 -Wall -Wextra -Wpedantic \
    -I include -pthread src/main.cpp -o fx_pricer

./fx_pricer
```

### Run tests

```bash
g++-14 -std=c++23 -fcoroutines -O2 -Wall -Wextra -Wpedantic \
    -I include -pthread tests/test_pricer.cpp -o test_pricer
./test_pricer
# Expected: 92/92 tests passed
```

### CMake (alternative)

```bash
cmake -B build && cmake --build build
ctest --test-dir build

# Python bindings (optional, requires pybind11)
cmake -DBUILD_PYTHON=ON -B build && cmake --build build
```

### What you'll see

The pricer runs 19 demos in order:

| Section | What it shows |
|---------|---------------|
| FX Option Pricing | GK analytical price + all Greeks |
| Variance Reduction | Sobol, antithetic, control variate vs plain MC |
| Convergence | MC converging to BS as paths increase |
| Greeks | Pathwise delta, likelihood-ratio gamma |
| Threading | jthread pool speedup |
| Early Stop | Coroutine stopping at 90K paths vs 5M budget (816x faster) |
| Error Handling | `std::expected` error propagation |
| Implied Vol | Newton-Raphson IV solver round-trip |
| Arena | Arena allocator zero-malloc demo |
| Heston | Volatility smile across strikes |
| Heston MC | QE discretisation vs semi-closed-form |
| Barrier | Down-and-out call at different barrier levels |
| FD PDE | European vs American put, BS cross-validation |
| Algo Comparison | TWAP vs VWAP vs IS on 10M EURUSD |
| Urgency Sweep | IS with different λ values |
| TCA Report | Full post-trade cost breakdown |
| Alpha Signals | Spread + momentum signal scores over ticks |
| Performance | Sharpe, Sortino, Calmar, max drawdown |
| Machine Epsilon | IEEE 754 double precision boundary demo |
Here’s a README‑friendly markdown version that will render cleanly on GitHub while keeping the terminal flavour.

# FULL output of the pricer + engine:

C++20 FX Option Pricer & Execution Engine

```text
╔════════════════════════════════════════════════════════════════╗
║  C++20 FX Option Pricer & Execution Engine                     ║
║  Pricing: GK · MC · Heston · Barrier · FD PDE                  ║
║  Execution: TWAP · VWAP · IS · TCA · Alpha Signals             ║
╚════════════════════════════════════════════════════════════════╝
```

---

## EURUSD 3M 1.0900 — Garman–Kohlhagen vs Monte Carlo

```text
S=1.0850 K=1.0900 T=0.25 σ=7.50% r_d=4.35% r_f=2.50%

BS Call: [GK-Analytic] px=0.016144 se=0.00e+00 10µs
Δ=+0.504606 Γ=9.742094 V=0.0022 Θ=-0.0001 ρ=0.0013

BS Put:  [GK-Analytic] px=0.016115 se=0.00e+00 0µs
Δ=-0.489163 Γ=9.742094 V=0.0022 Θ=-0.0001 ρ=-0.0014

Parity: C-P=0.00002943  S·df_f-K·df_d=0.00002943  err=1.11e-16

MC Call: [MC-CV+AT 500K paths] px=0.016144 se=3.26e-05 17818µs
Δ=+0.504599 Γ=9.732620 V=0.0000 Θ=0.0000 ρ=0.0000

MC Put:  [MC-CV+AT 500K paths] px=0.016114 se=3.41e-05 17530µs
Δ=-0.489170 Γ=9.734035 V=0.0000 Θ=0.0000 ρ=0.0000

BS-MC diff: 6.13e-07 (0.0 SE)
```

---

## Variance Reduction

```text
Analytic: 10.450584

Method                              Price    Std Err    |Error|         µs
------------------------------------------------------------------------
Raw (Sobol only)                10.450556   0.032912   0.000027       6553
+ Antithetic                    10.450002   0.032907   0.000582       5477
+ Control Variate               10.450477   0.019358   0.000107       8406
+ Antithetic + CV (full)        10.450365   0.019358   0.000219       7260
```

---

## Convergence: Sobol QMC vs Pseudo‑Random

```text
Analytic: 10.450584

     Paths     MC Price      Std Err      |Error|         µs
------------------------------------------------------------
      1000    10.414230     0.272009     0.036353         32
      5000    10.442332     0.122176     0.008252        155
     10000    10.447996     0.086499     0.002588        486
     50000    10.450025     0.038710     0.000559       1593
    100000    10.450266     0.027374     0.000318       4127
    500000    10.450383     0.012243     0.000201      18493
   1000000    10.450481     0.008658     0.000103      39758
```

---

## Pathwise Greeks vs Closed‑Form

```text
     Greek     Analytic  MC Pathwise        Error
----------------------------------------------------
     Delta     0.636831     0.636825     0.000006
     Gamma     0.018762     0.018746     0.000016
```

---

## Threading Speedup (1M paths)

```text
 Threads        Price           µs      Speedup
------------------------------------------------
       1    10.450284        29020         1.00x
       2    10.450476        12949         2.24x
       4    10.450500         8881         3.27x
```

---

## Coroutine Early Termination (target SE = 0.05)

```text
With early stop:
  [MC-+AT+ES 90K paths] px=10.445628 se=4.90e-02 2784µs
  Δ=+0.636838 Γ=0.018720 V=0.0000 Θ=0.0000 ρ=0.0000

Without early stop:
  [MC-+AT 5000K paths] px=10.446172 se=6.58e-03 1944065µs
  Δ=+0.636832 Γ=0.018737 V=0.0000 Θ=0.0000 ρ=0.0000

Speedup: 698.3x (stopped early once SE < 0.05)
```

---

## Error Handling (std::expected)

```text
Neg spot:   [E1001] S<=0 (include/pricer/models/../engines/black_scholes.hpp:118)
Zero vol:   [E1003] σ<=0 (include/pricer/models/../engines/black_scholes.hpp:120)
Chain OK:   Validated: 10.4506
```

---

## Implied Vol Solver (Newton–Raphson)

```text
Target (from σ=7.50%): 0.01614438
Recovered: 7.500000% (err=4.16e-17)
```

---

## Arena Allocator

```text
Capacity: 1048576 bytes
Alloc 1000 doubles: 8000 bytes used, span size=1000
Alloc 5000 doubles: 48000 bytes used, span size=5000
s1 = 999
After reset: 0 bytes used
```

---

## Heston Stochastic Vol — Smile Generation

```text
v0=0.04 κ=1.5 θ=0.04 ξ=0.3 ρ=-0.7
Feller condition: satisfied

  Strike  Moneyness   Call Price     Impl Vol
----------------------------------------------
   80.00       0.80      22.4553       24.03%
   85.00       0.85      18.0055       22.96%
   90.00       0.90      13.8369       21.89%
   95.00       0.95      10.0677       20.83%
  100.00       1.00       6.8257       19.77%
  105.00       1.05       4.2259       18.73%
  110.00       1.10       2.3341       17.74%
  115.00       1.15       1.1260       16.85%
  120.00       1.20       0.4705       16.11%

Negative ρ → downside skew (OTM put IV > ATM > OTM call)
```

---

## Heston MC (QE Scheme) vs Semi‑Closed‑Form

```text
Semi-closed-form: 10.361869

     Paths    Steps     MC Price         SE    |Error|
------------------------------------------------------
     50000       50    10.410026     0.0537     0.0482
    100000       50    10.372209     0.0379     0.0103
    200000       50    10.388740     0.0268     0.0269
    500000       50    10.362968     0.0169     0.0011
```

---

## Barrier Options — Down‑and‑Out Call (Brownian Bridge)

```text
Vanilla call: 10.4506

 Barrier   DO Price         SE    % Vanilla
--------------------------------------------
    60.0    10.4679     0.0330       100.2%
    70.0    10.4464     0.0330       100.0%
    80.0    10.3303     0.0329        98.8%
    85.0     9.9635     0.0331        95.3%
    90.0     8.6510     0.0326        82.8%
    95.0     5.6430     0.0292        54.0%

Barrier closer to spot → more knock-out probability → cheaper
```

---

## FD PDE Solver — European & American Put

```text
              Method      Price      Delta      Gamma
-----------------------------------------------------
       BS European     5.5735    -0.3632   0.018762
  FD European (CN)     5.5727    -0.3618   0.018725
  FD American (CN)     6.0874    -0.4092   0.022923

FD vs BS error:         0.000831 (0.015%)
Early exercise premium: 0.514657

Grid refinement study (European call):
   N×M   FD Price      |Error|         µs
------------------------------------------
 50×50   10.417494     0.033089         54
100×100  10.448316     0.002268        201
200×200  10.449494     0.001089        827
400×400  10.449893     0.000691       3186
```

---

## Algo Execution: TWAP vs VWAP vs IS (10M EURUSD buy)

```text
Order: BUY 10000000 EURUSD over 60 min
Impact model: η=0.20 σ_daily=0.700%

Algorithm           Fills  Avg Price  Arriv bps   VWAP bps     Impact
---------------------------------------------------------------------
TWAP                   19    1.08499      -0.10      +0.55      +0.55
VWAP                   20    1.08499      -0.10      +0.46      +0.46
IS-Almgren             19    1.08504      +0.38      +1.01      +1.01
```

---

## IS Urgency Sweep: Impact vs Drift Trade‑off

```text
Urgency         Arriv bps   Impact bps   Timing bps
--------------------------------------------------
λ=0.1               +0.29        +0.92        -2.19
λ=0.5               +0.34        +0.97        -2.19
λ=1.0               +0.38        +1.01        -2.20
λ=2.0               +0.46        +1.08        -2.20
λ=5.0               +0.58        +1.20        -2.22

Low λ → front-loads (high impact, low drift risk)
High λ → spreads evenly (low impact, high drift risk)
```

---

## Full TCA Report — VWAP Algo

```text
[VWAP] TCA Report
Fills: 20 | Qty: 25000000 | Duration: 3420000ms | Fill rate: 100.0%
Avg fill: 1.08512 | Arrival mid: 1.08500 | Terminal mid: 1.08502

── Cost Breakdown (bps) ──
VWAP slippage:      +0.49
Spread cost:        +0.57
Impact cost:        +0.49

── Risk ──
Timing risk:        +0.17 bps
Participation:      0.5%
```

---

## Execution Signals: Spread + Momentum

```text
  Tick        Mid     Spread    SpreadSig  MomentumSig
------------------------------------------------------
     0    1.08500       0.99      +0.0000      +0.0000
    20    1.08499       0.80      +0.1303      -0.0008
    40    1.08500       1.06      -0.1357      +0.0004
    60    1.08501       0.83      +0.1171      +0.0010
    80    1.08501       1.05      -0.1148      +0.0011
   100    1.08501       0.80      +0.1324      -0.0000
   120    1.08502       1.03      -0.1001      +0.0012
   140    1.08503       0.87      +0.0949      +0.0011
   160    1.08503       0.90      +0.0591      +0.0006
   180    1.08502       1.06      -0.1122      -0.0006
```

---

## Performance Analytics

```text
--- Synthetic strategy (252 days, μ=8bp/day) ---
Performance Report (252 periods)
Annualized return:      +0.6038  (+60.38%)
Annualized vol:          0.2062  (20.62%)
Sharpe ratio:           +2.3951
Sortino ratio:          +3.9347
Max drawdown:            0.0852  (8.52%)
Max DD duration:             47  periods
Calmar ratio:           +7.0837

---
```
## Demo Status

```text
All demos complete.
```
---


