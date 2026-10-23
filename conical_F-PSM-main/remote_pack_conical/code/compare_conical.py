"""compare_conical -- score the Lumerical conical tables against the conical referee and the two
closed forms, and draw one figure per case.

    python compare_conical.py                 # every results/lumerical_conical_<case>.csv present
    python compare_conical.py lumerical_R     # one case

Reads  results/reference_conical.csv          (benchmark PC: referee + closed forms, cases_conical.py)
       results/lumerical_conical_<case>.csv   (work PC: lumerical_rcwa_conical.py)
Writes results/summary_conical.csv     one row per (case, phi): rms and max of |eta_lum - eta_ref|
                                       for s and p incidence, the closed forms' errors against the
                                       referee on the same rows, power-sum extremes, the amplitude
                                       convention found, and -- when it is found -- the
                                       polarisation-uniform errors (report Eq. 65) of Lumerical
                                       against the referee and of the corrected closed form
                                       against Lumerical
       figures/conical_<case>.png      referee (line), corrected (dotted), Lumerical (markers)

Amplitude convention (UNTESTED semantics of grating_characterization).  On every row two readings
are tried for the signal order's complex amplitudes: "in-first" (Tsp = s in, p out) and "out-first"
(Tsp = p in, s out), each as a power-normalised set (|amp|^2 summing to the order's recorded power)
or a field-normalised one (the same after one common factor).  The reading that reproduces the
recorded s and p powers to 1e-5 on every row is reported; otherwise "unresolved" and the Jones
columns stay empty.  The sign of Lumerical's p-polarised source relative to the referee's e_p is
a further convention the powers cannot fix; both signs are tried and the smaller Jones error is
reported with the sign used (ASSUMED until a phase-sensitive check exists).
"""
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
AMP = ("ss", "sp", "ps", "pp")
deg = np.pi / 180.0


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def key(r):
    return (r["case"], round(float(r["phi_deg"]), 4), round(float(r["theta_deg"]), 4))


def herm_err(A, B):
    D = A - B
    D11, D22, D12 = D[0, 0].real, D[1, 1].real, D[0, 1]
    return float(abs(0.5 * (D11 + D22)) + np.sqrt((0.5 * (D11 - D22)) ** 2 + abs(D12) ** 2))


def amplitude_convention(rows, tr, ls, lp):
    """Returns (label, factor, in_first) or (\"unresolved\", None, None)."""
    a = {}
    for k in AMP:
        a[k] = np.array([complex(f(r.get(f"{tr}{k}_re")), f(r.get(f"{tr}{k}_im"))) for r in rows])
        if not np.isfinite(a[k]).all():
            return "unresolved (no amplitudes recorded)", None, None
    for in_first, (cs, cp) in (("in-first", (("ss", "sp"), ("ps", "pp"))), ("out-first", (("ss", "ps"), ("sp", "pp")))):
        ps = np.abs(a[cs[0]]) ** 2 + np.abs(a[cs[1]]) ** 2      # power for s incidence under this reading
        pp = np.abs(a[cp[0]]) ** 2 + np.abs(a[cp[1]]) ** 2
        if np.max(np.abs(ps - ls)) < 1e-5 and np.max(np.abs(pp - lp)) < 1e-5:
            return f"{in_first}, power-normalised", 1.0, in_first == "in-first"
        ok = (ps > 1e-6) & (pp > 1e-6)
        if ok.sum() >= 3:
            ratio = np.concatenate([ls[ok] / ps[ok], lp[ok] / pp[ok]])
            if np.std(ratio) < 1e-5 * max(1.0, abs(np.mean(ratio))):
                return f"{in_first}, field-normalised (x {np.mean(ratio):.6f})", float(np.mean(ratio)), in_first == "in-first"
    return "unresolved", None, None


def jones_from_row(r, tr, factor, in_first):
    g = {k: complex(f(r[f"{tr}{k}_re"]), f(r[f"{tr}{k}_im"])) for k in AMP}
    s = np.sqrt(factor)
    if in_first:
        J = np.array([[g["ss"], g["ps"]], [g["sp"], g["pp"]]])     # [out, in]
    else:
        J = np.array([[g["ss"], g["sp"]], [g["ps"], g["pp"]]])
    return s * J


def main(argv):
    import cases_conical as cc
    from rcwa_conical import solve as referee_solve
    from fpsm_conical import corrected
    geom = {}
    for c in cc.bench_cases() + cc.david_cases():
        geom[(c["case"], round(c["phi_deg"], 4))] = c
    ref = {key(r): r for r in read_csv(os.path.join(RES, "reference_conical.csv"))}
    names = argv or sorted({r["case"] for r in ref.values()})
    out_rows = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        plt = None
    for name in names:
        path = os.path.join(RES, f"lumerical_conical_{cc.file_stem(name)}.csv")
        if not os.path.isfile(path):
            print(f"{name}: no Lumerical table yet ({os.path.basename(path)})")
            continue
        lum = read_csv(path)
        phis = sorted({round(float(r["phi_deg"]), 4) for r in lum})
        if plt:
            fig, axes = plt.subplots(1, len(phis), figsize=(4.2 * len(phis), 3.6), squeeze=False)
        for ax_i, phi in enumerate(phis):
            rows = [r for r in lum if round(float(r["phi_deg"]), 4) == phi]
            rows = [r for r in rows if key(r) in ref]
            if not rows:
                print(f"{name} phi {phi}: no matching reference rows")
                continue
            th = np.array([float(r["theta_deg"]) for r in rows])
            ls, lp = np.array([f(r["eta_s_lum"]) for r in rows]), np.array([f(r["eta_p_lum"]) for r in rows])
            rs = np.array([f(ref[key(r)]["eta_s"]) for r in rows])
            rp = np.array([f(ref[key(r)]["eta_p"]) for r in rows])
            cs = np.array([f(ref[key(r)]["eta_s_corr"]) for r in rows])
            cp = np.array([f(ref[key(r)]["eta_p_corr"]) for r in rows])
            ec = np.array([f(ref[key(r)]["err_corrected"]) for r in rows])
            el = np.array([f(ref[key(r)]["err_leading"]) for r in rows])
            ps = np.array([f(r["power_s_lum"]) for r in rows])
            pp = np.array([f(r["power_p_lum"]) for r in rows])
            rms_s, rms_p = np.sqrt(np.mean((ls - rs) ** 2)), np.sqrt(np.mean((lp - rp) ** 2))
            max_s, max_p = np.max(np.abs(ls - rs)), np.max(np.abs(lp - rp))
            c = geom[(name, phi)]
            tr = "T" if c["tau"] == +1 else "R"
            label, factor, in_first = amplitude_convention(rows, tr, ls, lp)
            e_lum_ref, e_corr_lum, psign = "", "", ""
            if factor is not None:
                best = None
                for sgn in (+1.0, -1.0):
                    worst_lr, worst_cl = 0.0, 0.0
                    for r in rows:
                        Jl = jones_from_row(r, tr, factor, in_first) @ np.diag([1.0, sgn])
                        theta = float(r["theta_deg"]) * deg
                        R = referee_solve(c["k0"], c["n0"], c["n1"], c["K"], c["Psi"], theta, phi * deg, c["d"], M=10)
                        Jr = R["J_T"] if c["tau"] == +1 else R["J_R"]
                        worst_lr = max(worst_lr, herm_err(Jl.conj().T @ Jl, Jr.conj().T @ Jr))
                        if ref[key(r)]["eta_s_corr"] not in ("", None) and not str(ref[key(r)]["err_corrected"]).startswith("refused"):
                            C = corrected(c["k0"], c["n0"], c["n1"], c["K"], c["Psi"], theta, phi * deg, c["d"], c["tau"])
                            worst_cl = max(worst_cl, herm_err(C.conj().T @ C, Jl.conj().T @ Jl))
                    if best is None or worst_lr < best[0]:
                        best = (worst_lr, worst_cl, sgn)
                e_lum_ref, e_corr_lum, psign = best
            out_rows.append(dict(case=name, phi_deg=phi, n=len(rows), rms_s=rms_s, max_s=max_s, rms_p=rms_p, max_p=max_p,
                                 corr_max_err_vs_referee=np.nanmax(ec), lead_max_err_vs_referee=np.nanmax(el),
                                 power_s_min=ps.min(), power_s_max=ps.max(), power_p_min=pp.min(), power_p_max=pp.max(),
                                 amplitude_convention=label, p_source_sign=psign,
                                 jones_err_lumerical_vs_referee=e_lum_ref, jones_err_corrected_vs_lumerical=e_corr_lum))
            jtxt = f"   Jones: Lumerical vs referee {e_lum_ref:.2e}, corrected vs Lumerical {e_corr_lum:.2e} (p sign {psign:+.0f})" if factor else ""
            print(f"{name:14s} phi {phi:6.1f}: n {len(rows):3d}  s rms {rms_s:.2e} max {max_s:.2e}   p rms {rms_p:.2e} max {max_p:.2e}   "
                  f"corrected max err {np.nanmax(ec):.2e}  leading {np.nanmax(el):.2e}   power {min(ps.min(), pp.min()):.6f}..{max(ps.max(), pp.max()):.6f}   "
                  f"amp: {label}{jtxt}")
            if plt:
                ax = axes[0, ax_i]
                ax.plot(th, rs, "-", color="C0", label="referee s")
                ax.plot(th, rp, "-", color="C3", label="referee p")
                ax.plot(th, cs, ":", color="C0", label="corrected s")
                ax.plot(th, cp, ":", color="C3", label="corrected p")
                ax.plot(th, ls, "^", color="C0", ms=4, mfc="none", label="Lumerical s")
                ax.plot(th, lp, "v", color="C3", ms=4, mfc="none", label="Lumerical p")
                ax.set_title(f"{name}, phi = {phi:.0f} deg")
                ax.set_xlabel("theta (deg)")
                ax.set_ylabel("order-1 efficiency")
                if ax_i == 0:
                    ax.legend(fontsize=7)
        if plt:
            fig.tight_layout()
            os.makedirs(FIG, exist_ok=True)
            fig.savefig(os.path.join(FIG, f"conical_{cc.file_stem(name)}.png"), dpi=140)
            plt.close(fig)
    if out_rows:
        with open(os.path.join(RES, "summary_conical.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
            w.writeheader()
            w.writerows(out_rows)
        print(f"-> results/summary_conical.csv ({len(out_rows)} rows)")


if __name__ == "__main__":
    main(sys.argv[1:])
