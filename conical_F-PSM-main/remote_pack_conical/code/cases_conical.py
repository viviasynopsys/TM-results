"""cases_conical -- the conical replay cases for the Lumerical leg.

Two sets, one definition each, so the referee tables, the closed forms and the Lumerical driver
can never drift apart:

  BENCH   the four benchmark gratings of ../../2026-10-16/code/cases.py that the Lumerical model
          already builds (lumerical_R, lumerical_T, zemax_kba_R, unslanted_T; the rot5 substitute
          is not needed, Lumerical handles K_z = 0 through the unslanted bypass), replayed on the
          conical Bragg cone: for each azimuth phi the polar angle theta_B(phi) that satisfies the
          report's Eq. (67) with the grating's own K and Psi, and an angle sweep theta_B +- 6 deg
          in 0.25 deg steps at the design wavelength.
  DAVID   the report's three Table 1 families at phi = 150 deg, n1 = 0.024 (mu = 0.016),
          d = 0.16/n1 = 6.667 um, u in {-2 .. 2}: the gratings the certificates were issued for,
          built in the Lumerical model from two recording beams that give the same K and Psi.

Conventions are the report's (basis.py): theta >= 0 from the surface normal, phi from +x' toward
+z, so a planar replay at negative benchmark theta is (|theta|, phi = 180 deg).  Lumerical's own
theta and phi follow its simulation axes (x tangential, z normal, phi from +x toward +y); the
driver fits the theta sign on the planar sweep as in September and applies the same map to phi.
Efficiencies are invariant under phi -> -phi (z-mirror), so only cos phi is physical.

    python cases_conical.py          # prints the cone table and writes ../results/reference_conical.csv
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "2026-10-16", "code")))
from cases import CASES          # noqa: E402

from basis import geometry        # noqa: E402
from rcwa_conical import solve, jones_error, eta_of   # noqa: E402
from fpsm_conical import leading, corrected, principal_open   # noqa: E402

deg = np.pi / 180.0
HALF_SWEEP, STEP = 6.0, 0.25
AZIMUTHS = (100.0, 120.0, 150.0, 180.0, 30.0, 45.0, 60.0, 90.0)   # tried in this order; kept where the cone exists


def bragg_theta_with(beta, K, Psi, phi):
    """Eq. (67) solved for theta at fixed K, Psi, phi: 2 beta [cos th cos Psi - sin th cos phi sin Psi] = K.
    Returns theta in (0, 80 deg) or None.  K and beta in one unit; beta = 2 pi n0 / lambda."""
    A, B = np.cos(Psi), np.cos(phi) * np.sin(Psi)
    R = np.hypot(A, B)
    x = K / (2 * beta * R)
    if abs(x) > 1.0:
        return None
    chi = np.arctan2(B, A)
    for th in (np.arccos(x) - chi, -np.arccos(x) - chi):
        if 0.0 < th < 80 * deg:
            return th
    return None


def bench_cases():
    out = []
    for name in ("lumerical_R", "lumerical_T", "zemax_kba_R", "unslanted_T"):
        c = CASES[name]
        k0, n0, n1, K, Psi, d = 2 * np.pi / c["lam_um"], c["n0"], c["n1"], c["K"] * 1e3, c["psi"], c["d_um"]
        beta = k0 * n0
        tau = -1 if name.endswith("_R") else +1
        kept = []
        for phi_deg in AZIMUTHS:
            th = bragg_theta_with(beta, K, Psi, phi_deg * deg)
            if th is None or th / deg < 3.0 or th / deg > 70.0:
                continue
            kept.append((phi_deg, th / deg))
            if len(kept) == 4:
                break
        kept.sort()
        for phi_deg, thb in kept:
            thetas = np.round(np.arange(thb - HALF_SWEEP, thb + HALF_SWEEP + 1e-9, STEP), 6)
            thetas = thetas[thetas > 0.5]
            out.append(dict(set="BENCH", case=name, tau=tau, k0=k0, n0=n0, n1=n1, K=K, Psi=Psi, d=d,
                            lam_um=c["lam_um"], phi_deg=phi_deg, theta_bragg_deg=thb, thetas_deg=thetas))
    return out


def david_cases():
    import sweep_david as sd
    out = []
    k0 = 2 * np.pi / sd.LAM
    n1 = 0.024
    for fam, (tau, Psi_deg, th0) in sd.FAMILIES.items():
        Psi, phi = Psi_deg * deg, 150.0 * deg
        K = sd.period_K(th0 * deg, phi, Psi)
        # the seven detunings u of the report's sweep plus a 0.3 deg grid over +-3 deg, for a curve
        thetas = np.unique(np.round(np.concatenate([[th0 + u * n1 / 0.02 for u in sd.US],
                                                    np.arange(th0 - 3.0, th0 + 3.0 + 1e-9, 0.3)]), 6))
        out.append(dict(set="DAVID", case=fam.replace(" ", "_"), tau=tau, k0=k0, n0=sd.N0, n1=n1, K=K, Psi=Psi,
                        d=0.16 / n1, lam_um=sd.LAM, phi_deg=150.0, theta_bragg_deg=th0, thetas_deg=thetas))
    return out


def file_stem(name):
    """File-name stem of a case.  Windows file names are case-insensitive, so the DAVID family
    'Unslanted_T' would overwrite the BENCH grating 'unslanted_T'; it gets a '_DAVID' stem.
    Added 2026-10-23 on the work PC (the first run of both overwrote the BENCH table)."""
    bench = {c["case"] for c in bench_cases()}
    return name + "_DAVID" if name not in bench and name.lower() in {b.lower() for b in bench} else name


def recording_beams(K, Psi, beta):
    """Two unit recording beams inside the medium, (x, z) = (tangential, normal) as the Lumerical
    model wants them (k1 = (sin a1, cos a1), k2 = (sin a2, +-cos a2)), whose difference beta (k1 - k2)
    is the grating vector (K_x, -K_y) in the report's convention (basis.py: order m carries
    (alpha + m K_x, sigma gamma_m) and the depth phase exp(i m K_y y), i.e. k_1 = k_0 + (K_x, -K_y)).
    Returns (angle1_deg, angle2_deg, holo_type)."""
    Kv = np.array([K * np.sin(Psi), -K * np.cos(Psi)]) / beta      # (x, z)
    if np.linalg.norm(Kv) > 2.0:
        raise ValueError("|K| > 2 beta: no recording beams at this wavelength")
    m = np.array([-Kv[1], Kv[0]])
    m = m / np.linalg.norm(m) * np.sqrt(1.0 - 0.25 * Kv @ Kv)
    if m[1] < 0:
        m = -m
    k1, k2 = 0.5 * Kv + m, -0.5 * Kv + m
    if k1[1] < 0:
        # cos(K.r) is even in K: swapping the beams flips K_vec and leaves the grating unchanged,
        # so beam 1 can always be the forward one (the model writes k1 = (sin a1, 0, +cos a1))
        k1, k2 = k2, k1
    a1 = np.degrees(np.arctan2(k1[0], k1[1]))
    if k2[1] < 0:
        return a1, np.degrees(np.arctan2(k2[0], -k2[1])), "R"
    return a1, np.degrees(np.arctan2(k2[0], k2[1])), "T"


def check_recording_beams(K, Psi, beta, a1, a2, typ, tol=1e-9):
    """The model's own construction: k1 = (sin a1, cos a1), k2 = (sin a2, +-cos a2); beta (k1 - k2)
    must be +-(K_x, -K_y).  Returns the relative residual."""
    k1 = np.array([np.sin(np.radians(a1)), np.cos(np.radians(a1))])
    k2 = np.array([np.sin(np.radians(a2)), (-1.0 if typ == "R" else 1.0) * np.cos(np.radians(a2))])
    Kv = beta * (k1 - k2)
    want = np.array([K * np.sin(Psi), -K * np.cos(Psi)])
    return min(np.linalg.norm(Kv - want), np.linalg.norm(Kv + want)) / K


def reference_rows(case, M=10):
    rows = []
    k0, n0, n1, K, Psi, d, tau = (case[k] for k in ("k0", "n0", "n1", "K", "Psi", "d", "tau"))
    phi = case["phi_deg"] * deg
    for th in case["thetas_deg"]:
        theta = th * deg
        R = solve(k0, n0, n1, K, Psi, theta, phi, d, M=M)
        J = R["J_T"] if tau == +1 else R["J_R"]
        G = geometry(k0, n0, K, Psi, theta, phi, [-1, 2])
        closed = bool(G["gsq"][0] < 0 and G["gsq"][1] < 0)
        row = dict(set=case["set"], case=case["case"], phi_deg=case["phi_deg"], theta_deg=th, lambda_um=case["lam_um"],
                   eta_s=eta_of(J, [1, 0]), eta_p=eta_of(J, [0, 1]),
                   eta_s_to_p=float(abs(J[1, 0]) ** 2), eta_p_to_s=float(abs(J[0, 1]) ** 2),
                   power_s=R["sum_s"], power_p=R["sum_p"], neighbours_closed=closed,
                   eta_s_corr="", eta_p_corr="", err_corrected="", err_leading="")
        if principal_open(k0, n0, K, Psi, theta, phi):
            try:
                C = corrected(k0, n0, n1, K, Psi, theta, phi, d, tau)
                L = leading(k0, n0, n1, K, Psi, theta, phi, d, tau)
                row.update(eta_s_corr=eta_of(C, [1, 0]), eta_p_corr=eta_of(C, [0, 1]),
                           err_corrected=jones_error(C, J), err_leading=jones_error(L, J))
            except ValueError as e:
                row["err_corrected"] = f"refused: {e}"
        rows.append(row)
    return rows


def main():
    allc = bench_cases() + david_cases()
    print("set    case              tau  phi(deg)  theta_B(deg)  Lambda(um)  Psi(deg)  n1     d(um)   sweep")
    for c in allc:
        print(f"{c['set']:6s} {c['case']:16s} {c['tau']:+d}   {c['phi_deg']:6.1f}   {c['theta_bragg_deg']:8.3f}    "
              f"{2 * np.pi / c['K']:.6f}   {c['Psi'] / deg:7.2f}  {c['n1']:.3f}  {c['d']:7.3f}  {c['thetas_deg'][0]:.2f}..{c['thetas_deg'][-1]:.2f} ({len(c['thetas_deg'])})")
    rows = []
    for c in allc:
        rows += reference_rows(c)
        print(f"   referee done: {c['set']} {c['case']} phi {c['phi_deg']}", flush=True)
    res = os.path.join(HERE, "..", "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "reference_conical.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows -> results/reference_conical.csv")
    print("\nrecording beams for the DAVID set (Lumerical model inputs; angles in the medium):")
    for c in david_cases():
        beta = c["k0"] * c["n0"]
        a1, a2, typ = recording_beams(c["K"], c["Psi"], beta)
        res = check_recording_beams(c["K"], c["Psi"], beta, a1, a2, typ)
        print(f"   {c['case']:14s} angle1 {a1:9.4f}  angle2 {a2:9.4f}  type {typ}  w_rec {c['lam_um']} um   "
              f"K-vector residual {res:.1e}")


if __name__ == "__main__":
    main()
