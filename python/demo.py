#!/usr/bin/env python3
"""
FX Option Pricer — Python demo.

Demonstrates all key features via the pybind11 bindings:
  1. BS vs MC pricing comparison
  2. Heston smile generation
  3. Barrier option pricing
  4. FD vs BS comparison (European + American)

Build the module first:
  cmake -DBUILD_PYTHON=ON -B build && cmake --build build
  export PYTHONPATH=build  # or wherever fx_pricer.*.so ends up
"""

import fx_pricer as fp


def demo_bs_vs_mc():
    print("=" * 70)
    print("  1. BS vs MC pricing comparison")
    print("=" * 70)

    mkt = fp.MarketSnap(S=1.085, sigma=0.075, r_d=0.0435, r_f=0.025)
    call = fp.Contract(K=1.09, T=0.25, type=fp.OptType.Call)
    put = fp.Contract(K=1.09, T=0.25, type=fp.OptType.Put)

    bs = fp.BSEngine()
    bs_call = bs.price(call, mkt)
    bs_put = bs.price(put, mkt)
    print(f"  BS Call: {bs_call.price:.6f}")
    print(f"  BS Put:  {bs_put.price:.6f}")

    cfg = fp.MCConfig()
    cfg.n_paths = 500_000
    cfg.antithetic = True
    cfg.control_variate = True
    cfg.n_threads = 1
    mc = fp.MCEngine(r_d=0.0435, sigma=0.075, config=cfg)
    mc_call = mc.price(call, mkt)
    print(f"  MC Call: {mc_call.price:.6f} ± {mc_call.std_err:.2e}")
    print(f"  Diff:    {abs(bs_call.price - mc_call.price):.2e}")
    print()


def demo_heston_smile():
    print("=" * 70)
    print("  2. Heston smile generation")
    print("=" * 70)

    params = fp.HestonParams()
    params.v0 = 0.04
    params.kappa = 1.5
    params.theta = 0.04
    params.xi = 0.3
    params.rho = -0.7
    print(f"  Feller condition: {'satisfied' if params.feller_satisfied() else 'VIOLATED'}")

    S, r_d, r_f, T = 100.0, 0.05, 0.0, 0.5
    strikes = [S * m for m in [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]]

    mkt = fp.MarketSnap(S=S, sigma=0.20, r_d=r_d, r_f=r_f)  # sigma unused by Heston

    cfg = fp.MCConfig()
    cfg.n_paths = 200_000
    cfg.steps = 50
    cfg.antithetic = True
    cfg.n_threads = 1
    hmc = fp.HestonMCEngine(r_d=r_d, params=params, config=cfg)

    print(f"\n  {'Strike':>8s}  {'Moneyness':>10s}  {'MC Price':>10s}  {'MC SE':>10s}")
    print("  " + "-" * 44)
    for K in strikes:
        c = fp.Contract(K=K, T=T, type=fp.OptType.Call)
        r = hmc.price(c, mkt)
        print(f"  {K:>8.2f}  {K/S:>10.2f}  {r.price:>10.4f}  {r.std_err:>10.4f}")
    print()


def demo_barrier():
    print("=" * 70)
    print("  3. Barrier option pricing — Down-and-Out Call")
    print("=" * 70)

    mkt = fp.MarketSnap(S=100.0, sigma=0.20, r_d=0.05, r_f=0.0)

    cfg = fp.MCConfig()
    cfg.n_paths = 200_000
    cfg.steps = 100
    cfg.antithetic = True
    cfg.n_threads = 1
    bmc = fp.BarrierMCEngine(r_d=0.05, sigma=0.20, config=cfg)

    # Vanilla for reference
    vanilla = fp.BSEngine().price(fp.Contract(K=100.0, T=1.0, type=fp.OptType.Call), mkt)
    print(f"  Vanilla call: {vanilla.price:.4f}")
    print()

    barriers = [70.0, 80.0, 85.0, 90.0, 95.0]
    print(f"  {'Barrier':>8s}  {'DO Price':>10s}  {'SE':>10s}  {'% of Vanilla':>14s}")
    print("  " + "-" * 48)
    for B in barriers:
        bc = fp.BarrierContract(
            K=100.0, T=1.0, type=fp.OptType.Call,
            barrier=B, barrier_type=fp.BarrierType.DownAndOut
        )
        r = bmc.price_barrier(bc, mkt)
        pct = 100 * r.price / vanilla.price if vanilla.price > 0 else 0
        print(f"  {B:>8.1f}  {r.price:>10.4f}  {r.std_err:>10.4f}  {pct:>13.1f}%")
    print()


def demo_fd():
    print("=" * 70)
    print("  4. FD vs BS — European and American put")
    print("=" * 70)

    mkt = fp.MarketSnap(S=100.0, sigma=0.20, r_d=0.05, r_f=0.0)
    bs = fp.BSEngine()

    # European put
    eu_put = fp.Contract(K=100.0, T=1.0, type=fp.OptType.Put)
    bs_r = bs.price(eu_put, mkt)

    cfg = fp.FDConfig()
    cfg.n_space = 300
    cfg.n_time = 300
    cfg.theta = 0.5
    fd = fp.FDEngine(config=cfg)
    fd_r = fd.price(eu_put, mkt)

    print(f"  European Put BS:  {bs_r.price:.6f}  Δ={bs_r.greeks.delta:+.4f}")
    print(f"  European Put FD:  {fd_r.price:.6f}  Δ={fd_r.greeks.delta:+.4f}")
    print(f"  FD error:         {abs(bs_r.price - fd_r.price):.6f} ({100*abs(bs_r.price - fd_r.price)/bs_r.price:.3f}%)")
    print()

    # American put
    am_put = fp.Contract(K=100.0, T=1.0, type=fp.OptType.Put,
                          exercise=fp.Exercise.American)
    am_r = fd.price(am_put, mkt)
    premium = am_r.price - fd_r.price
    print(f"  American Put FD:  {am_r.price:.6f}  Δ={am_r.greeks.delta:+.4f}")
    print(f"  Early exercise premium: {premium:.6f}")
    print()


if __name__ == "__main__":
    print()
    print("  ╔════════════════════════════════════════════════════════════════╗")
    print("  ║  C++20 FX Option Pricer — Python Demo                        ║")
    print("  ╚════════════════════════════════════════════════════════════════╝")
    print()

    demo_bs_vs_mc()
    demo_heston_smile()
    demo_barrier()
    demo_fd()

    print("  All demos complete.")
    print()
