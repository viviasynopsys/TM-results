"""sweep_david -- reproduce the report's declared 504-point sweep (Sec. 5.A, Table 1-2) and its
Table 3 / Table 4 certificate points, scoring the transcribed closed forms (fpsm_conical.py)
against the independent conical referee (rcwa_conical.py) with the report's polarisation-uniform
error, Eq. (63)-(65).

    python sweep_david.py            # writes ../results/sweep_504.csv, prints Table 2 and the checks

Report values typed in below are the printed ones (Table 2, Table 3, Table 4, Sec. 5.C); they are
the oracle this script compares against.  lambda = 0.532 um, n0 = 1.5, phi_g = 0, d = 0.16/n1 um,
theta = theta0 + u (n1/0.02) deg, period fixed by Eq. (67) at (theta0, phi) and not retuned.
"""
import csv
import os
import sys
import time

import numpy as np

from basis import geometry
from rcwa_conical import solve, jones_error
from fpsm_conical import leading, corrected

deg = np.pi / 180.0
LAM, N0 = 0.532, 1.5
FAMILIES = {                       # name: (tau, Psi deg, theta0 deg)   -- Table 1
    "Slanted T": (+1, 80.0, 25.0),
    "Slanted R": (-1, 30.0, 20.0),
    "Unslanted T": (+1, 90.0, 25.0),
}
PHIS = (120.0, 135.0, 150.0, 165.0)
N1S = (0.003, 0.006, 0.012, 0.024, 0.048, 0.096)
US = (-2.0, -1.0, -0.35, 0.0, 0.35, 1.0, 2.0)
M_REF = 9                          # 19 retained orders, the report's comparison reference

TABLE2 = {                         # (family, model): maxima at mu = 0.002, 0.016, 0.064
    ("Slanted T", "L"): (2.205e-4, 1.990e-3, 1.407e-2),
    ("Slanted T", "C"): (1.966e-7, 4.029e-5, 6.289e-4),
    ("Slanted R", "L"): (1.032e-4, 8.640e-4, 6.365e-3),
    ("Slanted R", "C"): (2.714e-9, 1.658e-6, 1.378e-3),
    ("Unslanted T", "L"): (1.140e-4, 5.308e-3, 1.066e-1),
    ("Unslanted T", "C"): (4.771e-7, 2.745e-4, 1.527e-2),
}
TABLE3 = {                         # E_19 at phi = 150, u = 0, d = 0.16/n1; n1 = .006 .012 .024 .048
    "Slanted T": (2.287e-8, 1.669e-7, 1.366e-6, 9.774e-6),
    "Slanted R": (1.527e-8, 7.080e-8, 8.393e-7, 5.474e-6),
    "Unslanted T": (3.838e-8, 1.239e-7, 8.923e-7, 6.541e-6),
}
TABLE4 = {                         # E_inf at d = 2 um, phi = 150, u = 0; n1 = .006 .012 .024
    "Slanted T": (2.192e-9, 3.424e-8, 5.116e-7),
    "Slanted R": (1.047e-9, 1.651e-8, 2.501e-7),
    "Unslanted T": (2.396e-9, 3.870e-8, 6.392e-7),
}


def period_K(theta0, phi, Psi):
    """Eq. (67): bare first-order Bragg matching at the nominal direction."""
    beta = 2 * np.pi * N0 / LAM
    return 2 * beta * (np.cos(theta0) * np.cos(Psi) - np.sin(theta0) * np.cos(phi) * np.sin(Psi))


def point(fam, phi_deg, n1, u, d=None, M=M_REF):
    tau, Psi_deg, th0 = FAMILIES[fam]
    Psi, phi = Psi_deg * deg, phi_deg * deg
    k0 = 2 * np.pi / LAM
    K = period_K(th0 * deg, phi, Psi)
    theta = (th0 + u * n1 / 0.02) * deg
    d = 0.16 / n1 if d is None else d
    R = solve(k0, N0, n1, K, Psi, theta, phi, d, M=M)
    J = R["J_T"] if tau == +1 else R["J_R"]
    L = leading(k0, N0, n1, K, Psi, theta, phi, d, tau)
    C = corrected(k0, N0, n1, K, Psi, theta, phi, d, tau)
    eL = jones_error(L, J)
    eC = jones_error(C, J)
    G = geometry(k0, N0, K, Psi, theta, phi, [-1, 2])
    closed = bool(G["gsq"][0] < 0 and G["gsq"][1] < 0)
    return dict(family=fam, phi=phi_deg, n1=n1, mu=n1 / N0, u=u, theta_deg=theta / deg, d_um=d,
                Lambda_um=2 * np.pi / K, neighbours_closed=closed, err_leading=eL, err_corrected=eC,
                eta_s_rcwa=float(np.vdot(J[:, 0], J[:, 0]).real), eta_p_rcwa=float(np.vdot(J[:, 1], J[:, 1]).real),
                eta_s_corr=float(np.vdot(C[:, 0], C[:, 0]).real), eta_p_corr=float(np.vdot(C[:, 1], C[:, 1]).real),
                eta_s_lead=float(np.vdot(L[:, 0], L[:, 0]).real), eta_p_lead=float(np.vdot(L[:, 1], L[:, 1]).real),
                power_s=R["sum_s"], power_p=R["sum_p"])


def jones_error_op(CtC, J):
    """Eq. (65) with the closed-form side given as its efficiency operator C^H C."""
    D = CtC - J.conj().T @ J
    D11, D22, D12 = D[0, 0].real, D[1, 1].real, D[0, 1]
    return abs(0.5 * (D11 + D22)) + np.sqrt((0.5 * (D11 - D22)) ** 2 + abs(D12) ** 2)


def fmt(x):
    return f"{x:.3e}"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    res = os.path.join(here, "..", "results")
    os.makedirs(res, exist_ok=True)
    t0 = time.time()
    rows = []
    for fam in FAMILIES:
        for phi in PHIS:
            for n1 in N1S:
                for u in US:
                    rows.append(point(fam, phi, n1, u))
        print(f"  {fam}: done ({time.time() - t0:.0f} s)", flush=True)
    keys = list(rows[0])
    with open(os.path.join(res, "sweep_504.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} points -> results/sweep_504.csv")

    print("\nTable 2 reproduction: max over 28 (phi, u) points per cell.  printed | here | ratio")
    worst_ratio = 0.0
    for fam in FAMILIES:
        for model, key in (("L", "err_leading"), ("C", "err_corrected")):
            vals = []
            for n1 in (0.003, 0.024, 0.096):
                vals.append(max(r[key] for r in rows if r["family"] == fam and r["n1"] == n1))
            printed = TABLE2[(fam, model)]
            ratios = [v / p for v, p in zip(vals, printed)]
            worst_ratio = max(worst_ratio, max(abs(np.log(x)) for x in ratios))
            print(f"  {fam:12s} {model}  " + "   ".join(f"{fmt(p)} | {fmt(v)} | {x:.4f}" for p, v, x in zip(printed, vals, ratios)))
    print(f"  worst |ln ratio| over the 18 cells: {worst_ratio:.2e}")

    nclosed = sum(r["neighbours_closed"] for r in rows)
    print(f"\nneighbours -1, 2 both closed at {nclosed} points, open at {len(rows) - nclosed} (report: 401 / 103)")
    wc = max(rows, key=lambda r: r["err_corrected"])
    wl = max(rows, key=lambda r: r["err_leading"])
    print(f"worst corrected error {wc['err_corrected']:.4e} at {wc['family']}, phi {wc['phi']}, n1 {wc['n1']}, u {wc['u']}"
          f"   (report: 1.5264e-2 at Unslanted T, phi 120, n1 0.096, u 2)")
    print(f"worst leading error   {wl['err_leading']:.4e} at {wl['family']}, phi {wl['phi']}, n1 {wl['n1']}, u {wl['u']}"
          f"   (report: 1.0651e-1)")
    sub = [r for r in rows if r["mu"] <= 0.016 + 1e-12 and r["neighbours_closed"]]
    sub_open = [r for r in rows if r["mu"] <= 0.016 + 1e-12]
    print(f"subset mu <= 0.016 and neighbours closed: {len(sub)} points, max corrected {max(r['err_corrected'] for r in sub):.4e}"
          f" (report 268 points, 4.029e-5), max leading {max(r['err_leading'] for r in sub):.4e} (report 1.990e-3)")
    print(f"subset mu <= 0.016 incl. open neighbours: max corrected {max(r['err_corrected'] for r in sub_open):.4e} (report 2.745e-4)")
    print(f"power conservation: worst |1 - sum| = {max(max(abs(1 - r['power_s']), abs(1 - r['power_p'])) for r in rows):.2e}")

    print("\nTable 3: certified E_19 (upper bounds, 4 s.f.) against the observed error here, phi = 150, u = 0")
    print("  family        n1     E_19 printed   e(C, J_19) here   here/E_19")
    for fam in FAMILIES:
        for n1, E in zip((0.006, 0.012, 0.024, 0.048), TABLE3[fam]):
            r = point(fam, 150.0, n1, 0.0)
            print(f"  {fam:12s} {n1:.3f}   {E:.3e}      {r['err_corrected']:.4e}        {r['err_corrected'] / E:.4f}")

    print("\nTable 4: full-Maxwell E_inf (d = 2 um) against the observed error here with 27 orders (M = 13)")
    print("  family        n1     E_inf printed  e(C, J_27) here   here/E_inf")
    for fam in FAMILIES:
        for n1, E in zip((0.006, 0.012, 0.024), TABLE4[fam]):
            r = point(fam, 150.0, n1, 0.0, d=2.0, M=13)
            print(f"  {fam:12s} {n1:.3f}   {E:.3e}      {r['err_corrected']:.4e}        {r['err_corrected'] / E:.4f}")
    print(f"\ntotal {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
