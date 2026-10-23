"""test_conical -- the gates of the conical F-PSM test.  Run from code/:

    python test_conical.py            # prints every gate with its number, PASS/FAIL, and a summary
                                      # writes ../evidence/conical_gates.txt

G1  conical referee at phi = 0 reproduces the gated planar referee (../2026-10-16/code/rcwa.py) in
    both polarisations, every open order, R and T, on the five benchmark gratings (Bragg and +4 deg)
G2  power conservation of the conical referee at conical points (the five gratings x 5 azimuths,
    both polarisations) and over the report's 504-point sweep
G3  the report's Moharam-Gaylord Fig. 5 recomputation: E perp K first maximum at 61 um, H perp K at
    122 um, both above 0.99996 (report Sec. 5.B)
G4  Eq. (11) as printed (P_22 = -q A calA/k0) breaks power conservation; the derived order conserves it
G5  Table 2: the 18 cell maxima of the 504-point sweep, both models
G6  Table 3: the twelve finite-system certificates E_19 bound the observed e(C, J_19) and are tight
G7  Table 4: the nine full-Maxwell certificates bound the observed e(C, J_27)
G8  Sec. 5.C numbers: worst point, 401/103 split, 268-point subset maxima
G9  leading_scalar (Eq. 58-62) equals the four-state leading model at exact Bragg matching
G10 mirror symmetries of the referee: phi -> -phi and (phi -> pi - phi, Psi -> -Psi) leave efficiencies
G11 informational: the corrected model at phi = 0 on the five benchmark gratings -- the s channel
    against the September TE F-PSM (fpsm.py) and the referee, the p channel (an F-PSM TM result
    the 2026-10-16 work had no formula for) against the TM referee
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.normpath(os.path.join(HERE, "..", "..", "2026-10-16", "code"))
sys.path.insert(0, BENCH)
import rcwa as rcwa_planar                       # noqa: E402  (the gated 2026-10-16 referee)
import fpsm as fpsm_planar                       # noqa: E402
from cases import CASES                          # noqa: E402

from basis import geometry, e_vectors            # noqa: E402
from rcwa_conical import solve, jones_error, eta_of   # noqa: E402
from fpsm_conical import leading, leading_scalar, corrected   # noqa: E402
import sweep_david as sd                         # noqa: E402

deg = np.pi / 180.0
LOG = []
FAILS = []


def gate(name, ok, detail):
    line = f"{'PASS' if ok else 'FAIL'}  {name}: {detail}"
    LOG.append(line)
    print(line, flush=True)
    if not ok:
        FAILS.append(name)


def bench_geometry(c):
    """cases.py (nm) -> report conventions (um): K, Psi, d, k0."""
    return 2 * np.pi / c["lam_um"], c["n0"], c["n1"], c["K"] * 1e3, c["psi"], c["d_um"]


# ----------------------------------------------------------------------------- G1
def g1():
    worst = 0.0
    for name, c in CASES.items():
        k0, n0, n1, K, Psi, d = bench_geometry(c)
        for th in (c["theta_bragg_deg"], c["theta_bragg_deg"] + 4.0):
            # the report's theta is >= 0 with phi carrying the sign; cases.py's theta may be negative
            theta, phi = (th * deg, 0.0) if th >= 0 else (-th * deg, np.pi)
            rc = solve(k0, n0, n1, K, Psi, theta, phi, d, M=10)
            for pol, j in (("TE", 0), ("TM", 1)):
                p = rcwa_planar.solve(2 * np.pi / c["lam_nm"], n0, n1, c["K"], Psi, th * deg + Psi, c["d_nm"], L=10, pol=pol)
                worst = max(worst, np.abs(rc["DE_R"][:, j] - p["DE_R"][::-1]).max(),
                            np.abs(rc["DE_T"][:, j] - p["DE_T"][::-1]).max())
    gate("G1 planar reduction", worst <= 1e-10, f"max |dDE| over 5 gratings x 2 angles x 2 pols x all orders = {worst:.2e} (gate 1e-10)")


# ----------------------------------------------------------------------------- G2
def g2():
    worst = 0.0
    for name, c in CASES.items():
        k0, n0, n1, K, Psi, d = bench_geometry(c)
        for phi_deg in (30.0, 60.0, 90.0, 120.0, 150.0):
            for th in (10.0, 25.0, 45.0):
                rc = solve(k0, n0, n1, K, Psi, th * deg, phi_deg * deg, d, M=10)
                worst = max(worst, abs(1 - rc["sum_s"]), abs(1 - rc["sum_p"]))
    gate("G2a power, five gratings x 5 azimuths x 3 angles", worst <= 1e-10, f"max |1 - sum| = {worst:.2e} (gate 1e-10)")


# ----------------------------------------------------------------------------- G3
def g3():
    lam, Lam, Psi, phi = 0.5, 1.0 / 3.0, 90 * deg, 150 * deg
    theta = np.arcsin(1 / np.sqrt(3))
    n0 = 1.5
    k0, K = 2 * np.pi / lam, 2 * np.pi / Lam
    beta = k0 * n0
    eps_h = {0: 2.25, 1: 0.005, -1: 0.005}
    w = np.sqrt(2.25 ** 2 - 0.01 ** 2)
    r = 0.01 / (2.25 + w)
    inv = {h: (-r) ** abs(h) / w for h in range(-40, 41)}
    G = geometry(k0, n0, K, Psi, theta, phi, [0])
    am, q, gam = G["am"][0], G["q"], G["gam"][0]
    es, ep = e_vectors(am, q, gam, beta, +1.0)
    khat = np.array([am, gam, q]) / beta
    bE = np.array([-ep[0], es[0]])
    bE = bE / np.linalg.norm(bE)
    gs, gp = np.cross(khat, es), np.cross(khat, ep)
    bH = np.array([-gp[0], gs[0]])
    bH = bH / np.linalg.norm(bH)
    ds = np.arange(50.0, 130.01, 0.5)
    eE, eH = [], []
    for d in ds:
        J = solve(k0, n0, 0.0, K, Psi, theta, phi, d, M=7, eps_h=eps_h, inv_eps_h=inv)["J_T"]
        eE.append(eta_of(J, bE))
        eH.append(eta_of(J, bH))
    eE, eH = np.array(eE), np.array(eH)
    dE, vE = ds[np.argmax(eE[ds <= 90])], eE[ds <= 90].max()
    dH, vH = ds[np.argmax(eH)], eH.max()
    ok = abs(dE - 61.0) <= 1.0 and abs(dH - 122.0) <= 1.0 and vE > 0.99996 and vH > 0.99996
    gate("G3 Moharam-Gaylord Fig. 5 geometry", ok,
         f"E perp K first max at d = {dE:.1f} um, eta = {vE:.6f}; H perp K at {dH:.1f} um, eta = {vH:.6f} (report: 61, 122, > 0.99996)")


# ----------------------------------------------------------------------------- G4
def g4():
    c = CASES["lumerical_T"]
    k0, n0, n1, K, Psi, d = bench_geometry(c)
    worst_printed, worst_derived = 0.0, 0.0
    for phi_deg in (30.0, 90.0, 150.0):
        for order, acc in (("printed", "p"), ("derived", "d")):
            rc = solve(k0, n0, n1, K, Psi, 20 * deg, phi_deg * deg, d, M=10, p22_order=order)
            v = max(abs(1 - rc["sum_s"]), abs(1 - rc["sum_p"]))
            if order == "printed":
                worst_printed = max(worst_printed, v)
            else:
                worst_derived = max(worst_derived, v)
    gate("G4 Eq. (11) operator order", worst_printed > 1e-3 and worst_derived <= 1e-10,
         f"|1 - sum|: printed P_22 order {worst_printed:.2e}, derived order {worst_derived:.2e}")


# ----------------------------------------------------------------------------- G5-G8 (the sweep)
def g5_to_g8():
    t0 = time.time()
    rows = []
    for fam in sd.FAMILIES:
        for phi in sd.PHIS:
            for n1 in sd.N1S:
                for u in sd.US:
                    rows.append(sd.point(fam, phi, n1, u))
    print(f"   (504-point sweep: {time.time() - t0:.0f} s)", flush=True)
    worst_pow = max(max(abs(1 - r["power_s"]), abs(1 - r["power_p"])) for r in rows)
    gate("G2b power, 504-point sweep", worst_pow <= 1e-10, f"max |1 - sum| = {worst_pow:.2e} (gate 1e-10)")

    worst_rel, worst_cell = 0.0, ""
    for fam in sd.FAMILIES:
        for model, key in (("L", "err_leading"), ("C", "err_corrected")):
            for n1, printed in zip((0.003, 0.024, 0.096), sd.TABLE2[(fam, model)]):
                v = max(r[key] for r in rows if r["family"] == fam and r["n1"] == n1)
                rel = abs(v - printed) / printed
                if rel > worst_rel and abs(v - printed) > 1e-11:
                    worst_rel, worst_cell = rel, f"{fam} {model} n1={n1}: printed {printed:.3e}, here {v:.4e}"
    gate("G5 Table 2 (18 cells, both models)", worst_rel <= 2e-3,
         f"worst relative deviation {worst_rel:.2e} at {worst_cell or 'none'} (gate 2e-3 = 4-figure rounding + referee floor)")

    nclosed = sum(r["neighbours_closed"] for r in rows)
    wc = max(rows, key=lambda r: r["err_corrected"])
    wl = max(rows, key=lambda r: r["err_leading"])
    sub = [r for r in rows if r["mu"] <= 0.016 + 1e-12 and r["neighbours_closed"]]
    sub_open = [r for r in rows if r["mu"] <= 0.016 + 1e-12]
    mc, ml, mo = max(r["err_corrected"] for r in sub), max(r["err_leading"] for r in sub), max(r["err_corrected"] for r in sub_open)
    ok = (nclosed == 401 and len(sub) == 268
          and abs(wc["err_corrected"] - 1.5264e-2) / 1.5264e-2 <= 2e-3
          and wc["family"] == "Unslanted T" and wc["phi"] == 120.0 and wc["n1"] == 0.096 and wc["u"] == 2.0
          and abs(wl["err_leading"] - 1.0651e-1) / 1.0651e-1 <= 2e-3
          and abs(mc - 4.029e-5) / 4.029e-5 <= 2e-3 and abs(ml - 1.990e-3) / 1.990e-3 <= 2e-3
          and abs(mo - 2.745e-4) / 2.745e-4 <= 2e-3)
    gate("G8 Sec. 5.C numbers", ok,
         f"closed/open {nclosed}/{len(rows) - nclosed} (401/103); worst C {wc['err_corrected']:.4e} at {wc['family']} phi {wc['phi']:.0f} "
         f"n1 {wc['n1']} u {wc['u']} (1.5264e-2, Unslanted T 120 0.096 2); worst L {wl['err_leading']:.4e} (1.0651e-1); "
         f"268-subset C {mc:.4e} (4.029e-5), L {ml:.4e} (1.990e-3); mu<=0.016 incl. open C {mo:.4e} (2.745e-4)")

    bad = []
    for fam in sd.FAMILIES:
        for n1, E in zip((0.006, 0.012, 0.024, 0.048), sd.TABLE3[fam]):
            e = sd.point(fam, 150.0, n1, 0.0)["err_corrected"]
            if not (e <= E * (1 + 5e-4) + 1e-11 and e >= 0.99 * E):
                bad.append(f"{fam} n1={n1}: E_19 {E:.3e}, observed {e:.4e}")
    gate("G6 Table 3 certificates", not bad, "all twelve observed errors within [0.99, 1.0005] x E_19" if not bad else "; ".join(bad))

    bad = []
    for fam in sd.FAMILIES:
        for n1, E in zip((0.006, 0.012, 0.024), sd.TABLE4[fam]):
            e = sd.point(fam, 150.0, n1, 0.0, d=2.0, M=13)["err_corrected"]
            if not (e <= E * (1 + 5e-4) + 1e-11 and e >= 0.99 * E):
                bad.append(f"{fam} n1={n1}: E_inf {E:.3e}, observed {e:.4e}")
    gate("G7 Table 4 certificates", not bad, "all nine observed errors within [0.99, 1.0005] x E_inf" if not bad else "; ".join(bad))
    return rows


# ----------------------------------------------------------------------------- G9
def g9():
    worst = 0.0
    k0 = 2 * np.pi / sd.LAM
    for fam, (tau, Psi_deg, th0) in sd.FAMILIES.items():
        for phi_deg in sd.PHIS:
            for n1 in (0.006, 0.048):
                Psi, phi = Psi_deg * deg, phi_deg * deg
                K = sd.period_K(th0 * deg, phi, Psi)
                d = 0.16 / n1
                JL = leading(k0, sd.N0, n1, K, Psi, th0 * deg, phi, d, tau)
                S = leading_scalar(k0, sd.N0, n1, K, Psi, th0 * deg, phi, d, tau)
                worst = max(worst, sd.jones_error_op(S["CtC"], JL))
    gate("G9 scalar channels = four-state model at Bragg", worst <= 1e-12, f"max e(C_scalar, C_4state) at u = 0 over 24 points = {worst:.2e}")


# ----------------------------------------------------------------------------- G10
def g10():
    worst = 0.0
    for name in ("lumerical_T", "zemax_kba_R"):
        c = CASES[name]
        k0, n0, n1, K, Psi, d = bench_geometry(c)
        for phi_deg, th in ((40.0, 20.0), (130.0, 35.0)):
            a = solve(k0, n0, n1, K, Psi, th * deg, phi_deg * deg, d, M=8)
            b = solve(k0, n0, n1, K, Psi, th * deg, -phi_deg * deg, d, M=8)
            cc = solve(k0, n0, n1, K, -Psi, th * deg, (180 - phi_deg) * deg, d, M=8)
            for other in (b, cc):
                worst = max(worst, np.abs(a["DE_R"] - other["DE_R"]).max(), np.abs(a["DE_T"] - other["DE_T"]).max())
    gate("G10 mirror symmetries", worst <= 1e-10, f"max |dDE| under z-mirror and x'-mirror = {worst:.2e}")


# ----------------------------------------------------------------------------- G11 (informational)
def g11():
    LOG.append("INFO  G11 corrected conical model at phi = 0 on the five benchmark gratings (one length unit = um):")
    LOG.append("      case              theta     |TE: F-PSM(Sep) - RCWA|  |TE: conical C - RCWA|  |TE: conical C - F-PSM(Sep)|  |TM: conical C - RCWA|  TM RCWA")
    for name, c in CASES.items():
        k0, n0, n1, K, Psi, d = bench_geometry(c)
        tau = -1 if "_R" in name else +1
        for th in (c["theta_bragg_deg"], c["theta_bragg_deg"] + 2.0):
            theta, phi = (th * deg, 0.0) if th >= 0 else (-th * deg, np.pi)
            rc = solve(k0, n0, n1, K, Psi, theta, phi, d, M=10)
            J = rc["J_T"] if tau == +1 else rc["J_R"]
            try:
                C = corrected(k0, n0, n1, K, Psi, theta, phi, d, tau)
                te_c, tm_c = eta_of(C, [1, 0]), eta_of(C, [0, 1])
            except ValueError as e:
                LOG.append(f"      {name:16s} {th:8.3f}   corrected model refused: {e}")
                continue
            te_r, tm_r = eta_of(J, [1, 0]), eta_of(J, [0, 1])
            fp = fpsm_planar.efficiency(2 * np.pi / c["lam_nm"], n0, n1, c["K"], Psi, th * deg + Psi, c["d_nm"])
            fp = float(fp)
            LOG.append(f"      {name:16s} {th:8.3f}   {abs(fp - te_r):.2e}                 {abs(te_c - te_r):.2e}               "
                       f"{abs(te_c - fp):.2e}                     {abs(tm_c - tm_r):.2e}             {tm_r:.6f}")
    for line in LOG[-11:]:
        print(line)


def main():
    t0 = time.time()
    g1()
    g2()
    g3()
    g4()
    g5_to_g8()
    g9()
    g10()
    g11()
    summary = f"\n{len(FAILS)} failure(s) of 12 gates; {time.time() - t0:.0f} s; 2026-10-23"
    if FAILS:
        summary += ": " + ", ".join(FAILS)
    print(summary)
    LOG.append(summary)
    ev = os.path.join(HERE, "..", "evidence")
    os.makedirs(ev, exist_ok=True)
    with open(os.path.join(ev, "conical_gates.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
