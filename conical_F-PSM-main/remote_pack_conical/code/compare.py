"""F-PSM next to Kogelnik and RCWA, in the format of the Ansys articles' figure
(first-order efficiency vs replay angle), plus the wavelength sweep and the
error against the rigorous solution on a log axis.  TE or TM (2026-10-16).

    python compare.py [case ...]              TE: exactly the 30 September benchmark (regression)
    python compare.py --pol TM [case ...]     TM: referee, Kogelnik and the Ansys TM tables

For each case in ``cases.py`` (SUF = "" for TE, "_TM" for TM):
    ../results/compare_<case>SUF_angle.csv, _wavelength.csv     every curve, every point
    ../figures/compare_<case>SUF.png                            the side-by-side figure
and ../results/summarySUF.csv with rms / max errors against RCWA.

Overlaid when their CSVs exist in ../results (absence is reported, never papered over):
    zemax_<case>SUF.csv            OpticStudio native Kogelnik      (zemax_kogelnik_nsc.py)
    zemax_rcwa_<case>SUF.csv       OpticStudio RCWA hologram DLL    (zemax_rcwa_nsc.py)
    lumerical_rcwa_<case>SUF.csv   Lumerical RCWA solver            (lumerical_rcwa.py)
All three share one schema: sweep, source_tilt_x_deg (= -theta), lambda_um,
P_entry_side_W (reflected first order), P_far_side_W (transmitted first order).
The two Ansys RCWA tables are also scored: Python RCWA, Kogelnik and F-PSM against each.

Diagnostic (added 2026-09-18): the OpticStudio RCWA DLL expands the permittivity as
eps = n0^2 + 2 n0 dn cos(K.r) (no dn^2 terms), whereas the note / referee use n = n0 + n1 cos
exactly.  ``rcwa.signal_efficiency(..., eps_form="linear")`` reproduces the DLL; it is scored
and drawn as a thin dashed line so the reader can see how much of an Ansys deviation is that
definitional difference.  The referee (solid, black) is unchanged.

TM mode (2026-10-16).  F-PSM has no TM formulae (fpsm.TM_ABSENT), so its column is EMPTY, its
curve is absent and its scores are NaN -- nothing is substituted.  The domain flags written are
F-PSM's TE conditions evaluated on the same geometry (columns named fpsm_te_domain_*), kept so
the TM tables can be scored over the same point sets once a TM sheet exists.  Two extra columns:
eta_rcwa_L10_TE (the TE referee on the same points, for the polarisation contrast) and
rel_err_kogelnik = |eta_K - eta_ref| / max(eta_ref, 1e-3): TM efficiency collapses where the two
principal waves are orthogonal, and absolute error flatters every model there.  The TE path of
this script is byte-for-byte the 30 September one (regression: evidence/te_regression.txt).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fpsm
import kogelnik
import rcwa
from cases import CASES, deg

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES, FIG = os.path.join(ROOT, "results"), os.path.join(ROOT, "figures")
L_REF, L_CHECK = 10, 14
# PAPER mode (env COMPARE_PAPER=1): no suptitle (the caption carries it), vector SVG + PNG into
# figures/paper/, and the CSV record (compare_*.csv, summary.csv) is left untouched.  2026-09-30.
PAPER = os.environ.get("COMPARE_PAPER", "") not in ("", "0")
REL_FLOOR = 1e-3            # eta floor of the relative error (the benchmark's acceptance-gate scale)


def _pol(argv):
    if "--pol" in argv:
        i = argv.index("--pol")
        pol = argv[i + 1].upper()
        del argv[i:i + 2]
    else:
        pol = "TE"
    if pol not in rcwa.POLS:
        raise SystemExit(f"--pol must be TE or TM, not {pol!r}")
    return pol


POL = _pol(sys.argv)
TE = POL == "TE"
SUF = "" if TE else "_TM"

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C_KOG, C_FPSM, C_ZMX, C_ZRC, C_LUM = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"   # categorical slots 1-5
C_TEREF = "#9a9892"
ANSYS_RCWA = (("zemax_rcwa", "OpticStudio RCWA DLL", C_ZRC, "s"), ("lumerical_rcwa", "Lumerical RCWA", C_LUM, "^"))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb"})


def models(c, theta_deg, lam_nm):
    theta_deg, lam_nm = np.broadcast_arrays(np.asarray(theta_deg, float), np.asarray(lam_nm, float))
    k0 = 2 * np.pi / lam_nm
    # OpticStudio's convention (Kogelnik article, eq. 2): the traced "first order"
    # is the one for which K . k_playback > 0, so K flips sign where the replay
    # beam crosses the fringe planes.  K -> -K is psi -> psi + 180 deg here (same
    # grating, other order).  Only lumerical_T (theta > 14 deg) is affected.
    psi = np.where(np.cos(theta_deg * deg + c["psi"]) >= 0, c["psi"], c["psi"] + np.pi)
    thc = theta_deg * deg + psi
    a = (k0, c["n0"], c["n1"], c["K"], psi, thc, c["d_nm"])
    det = fpsm.fpsm_detail(*a)                                   # TE lossless: the geometry, flags and (TE) F-PSM
    eR, eT = rcwa.signal_efficiency(*a, L=L_REF, pol=POL)
    eR2, eT2 = rcwa.signal_efficiency(*a, L=L_CHECK, pol=POL)
    eRl, eTl = rcwa.signal_efficiency(*a, L=L_REF, eps_form="linear", pol=POL)   # the DLL's eps expansion
    refl = det["case"] == -1
    if not (np.all(refl[det["case"] != 0]) or not np.any(refl)):
        raise RuntimeError("geometry flag changes inside one sweep -- split the sweep")
    rc, rc2, rcl = (eR, eR2, eRl) if np.any(refl) else (eT, eT2, eTl)
    out = dict(rcwa=rc, rcwa_conv=float(np.max(np.abs(rc - rc2))), rcwa_epslin=rcl,
               in_domain=det["in_domain"], domain_fail=det["domain_fail"], applicable=det["case"] != 0,
               kogelnik=kogelnik.efficiency(*a, pol=POL), reflection=bool(np.any(refl)))
    if TE:
        out["fpsm"] = det["eta"]
    else:
        out["fpsm"] = np.full(rc.shape, np.nan)                  # fpsm.TM_ABSENT: no TM formulae, nothing substituted
        teR, teT = rcwa.signal_efficiency(*a, L=L_REF, pol="TE")
        out["rcwa_te"] = teR if out["reflection"] else teT
    return out


def load_zemax(name, reflection, prefix="zemax"):
    p = os.path.join(RES, f"{prefix}_{name}{SUF}.csv")
    if not os.path.isfile(p):
        return None
    rows = list(csv.DictReader(open(p)))
    col = "P_entry_side_W" if reflection else "P_far_side_W"
    ang = np.array([(-float(r["source_tilt_x_deg"]), float(r[col])) for r in rows if r["sweep"] == "angle"])
    wav = np.array([(float(r["lambda_um"]) * 1e3, float(r[col])) for r in rows if r["sweep"] == "wavelength"])
    return ang, wav


def valid(m):
    """Points scored: where F-PSM is defined (TE), i.e. the principal channel is open.  In TM the
    same geometric criterion (fpsm 'case' != 0) is used, so TE and TM are scored on identical sets."""
    return np.isfinite(m["fpsm"]) if TE else m["applicable"]


def stats(m, mask=None):
    mask = valid(m) if mask is None else mask
    out = {}
    for k in ("kogelnik", "fpsm"):
        e = (m[k] - m["rcwa"])[mask]
        out[k] = (float(np.sqrt(np.mean(e ** 2))), float(np.max(np.abs(e)))) if e.size and np.all(np.isfinite(e)) else (np.nan, np.nan)
    if not TE:
        r = (np.abs(m["kogelnik"] - m["rcwa"]) / np.maximum(m["rcwa"], REL_FLOOR))[mask]
        out["kogelnik_rel"] = (float(np.sqrt(np.mean(r ** 2))), float(np.max(r))) if r.size else (np.nan, np.nan)
    return out


def shade(ax, x, ok):
    """Grey band where F-PSM is OUTSIDE its stated working domain."""
    bad = ~ok
    if not bad.any():
        return
    edges = np.flatnonzero(np.diff(np.r_[0, bad.astype(int), 0]))
    for s, e in zip(edges[::2], edges[1::2]):
        ax.axvspan(x[s], x[min(e, len(x) - 1)], color="#000000", alpha=0.055, lw=0)


def write_table(path, xname, x, m):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        if TE:
            w.writerow([xname, "eta_rcwa_L10", "eta_kogelnik", "eta_fpsm", "fpsm_in_domain", "fpsm_domain_failmask",
                        "eta_rcwa_L10_epslinear"])
            for i in range(len(x)):
                w.writerow([f"{x[i]:.6f}", f"{m['rcwa'][i]:.12e}", f"{m['kogelnik'][i]:.12e}",
                            f"{m['fpsm'][i]:.12e}", int(m["in_domain"][i]), int(m["domain_fail"][i]),
                            f"{m['rcwa_epslin'][i]:.12e}"])
        else:
            w.writerow([xname, "eta_rcwa_L10", "eta_kogelnik", "eta_fpsm", "fpsm_te_domain_in", "fpsm_te_domain_failmask",
                        "eta_rcwa_L10_epslinear", "eta_rcwa_L10_TE", "rel_err_kogelnik", "pol"])
            for i in range(len(x)):
                rel = abs(m["kogelnik"][i] - m["rcwa"][i]) / max(m["rcwa"][i], REL_FLOOR)
                w.writerow([f"{x[i]:.6f}", f"{m['rcwa'][i]:.12e}", f"{m['kogelnik'][i]:.12e}",
                            "", int(m["in_domain"][i]), int(m["domain_fail"][i]),
                            f"{m['rcwa_epslin'][i]:.12e}", f"{m['rcwa_te'][i]:.12e}", f"{rel:.6e}", POL])


def one_case(name, c, summary):
    th, lam = c["theta_deg"], c["lam_sweep_nm"]
    A = models(c, th, c["lam_nm"])
    W = models(c, c["theta_bragg_deg"], lam)
    zm = load_zemax(name, A["reflection"])
    ansys = {pre: load_zemax(name, A["reflection"], pre) for pre, *_ in ANSYS_RCWA}
    for pre, tab in ansys.items():
        if tab is None:
            print(f"   [{name}] no {pre}_{name}{SUF}.csv -- that Ansys RCWA curve is absent from this figure")

    for tag, x, xname, m in (() if PAPER else (("angle", th, "theta_deg", A), ("wavelength", lam, "lambda_nm", W))):
        write_table(os.path.join(RES, f"compare_{name}{SUF}_{tag}.csv"), xname, x, m)

    fig, axs = plt.subplots(2, 2, figsize=(11, 6.6), sharex="col",
                            gridspec_kw=dict(height_ratios=[3, 2], hspace=0.08, wspace=0.16))
    order = "R" if A["reflection"] else "T"
    pol_letter = "s" if TE else "p"
    for j, (x, m, xl, z) in enumerate(((th, A, "Replay angle in the medium, θ [°]", None if zm is None else zm[0]),
                                       (lam, W, "Replay wavelength [nm]", None if zm is None else zm[1]))):
        ax, ex = axs[0, j], axs[1, j]
        for a_ in (ax, ex):
            a_.grid(True, color=GRID, lw=0.6)
            a_.set_axisbelow(True)
            shade(a_, x, m["in_domain"])
        if not TE:
            ax.plot(x, m["rcwa_te"], color=C_TEREF, lw=1.0, ls=(0, (2, 2)), label="RCWA, TE (same grating, for contrast)")
        ax.plot(x, m["rcwa"], color=INK, lw=2.0, label=f"RCWA (rigorous, L = {L_REF}){'' if TE else ', TM'}")
        ax.plot(x, m["kogelnik"], color=C_KOG, lw=1.6, label="Kogelnik" if TE else "Kogelnik, TM (κ cos 2θ_c)")
        if TE:
            ax.plot(x, m["fpsm"], color=C_FPSM, lw=1.6, ls=(0, (5, 2)), label="F-PSM")
        if z is not None and len(z):
            step = max(1, len(z) // 40)
            ax.plot(z[::step, 0], z[::step, 1], ls="none", marker="o", ms=4.2, mfc="none", mec=C_ZMX, mew=1.3,
                    label="OpticStudio native (Kogelnik)")
        for k_, (pre, lab, col, mk) in enumerate(ANSYS_RCWA):
            tab = None if ansys[pre] is None else ansys[pre][j]
            if tab is None or not len(tab):
                continue
            tab = tab[np.argsort(tab[:, 0])]
            step = max(1, len(tab) // 40)
            off = (k_ + 1) * step // 3            # stagger the marker sets so they do not sit on top of each other
            ax.plot(tab[off::step, 0], tab[off::step, 1], ls="none", marker=mk, ms=4.6, mfc="none", mec=col, mew=1.3, label=lab)
            ex.plot([], [], color=col, lw=1.2, ls=(0, (1, 1)), label=f"|{lab} − RCWA|")
            ex.semilogy(tab[:, 0], np.maximum(np.abs(tab[:, 1] - np.interp(tab[:, 0], x, m["rcwa"])), 1e-12),
                        color=col, lw=1.2, ls=(0, (1, 1)))
            if pre == "zemax_rcwa":              # the DLL's own eps expansion, for diagnosis only
                ex.plot([], [], color=col, lw=0.9, ls=(0, (4, 2)), alpha=0.8, label=f"|{lab} − RCWA with ε = n₀²+2n₀n₁·cos|")
                ex.semilogy(tab[:, 0], np.maximum(np.abs(tab[:, 1] - np.interp(tab[:, 0], x, m["rcwa_epslin"])), 1e-12),
                            color=col, lw=0.9, ls=(0, (4, 2)), alpha=0.8)
        ax.set_ylim(-0.02, 1.02)
        ax.set_ylabel(f"{order}{pol_letter} power, order +1")
        ex.semilogy(x, np.maximum(np.abs(m["kogelnik"] - m["rcwa"]), 1e-12), color=C_KOG, lw=1.6)
        if TE:
            ex.semilogy(x, np.maximum(np.abs(m["fpsm"] - m["rcwa"]), 1e-12), color=C_FPSM, lw=1.6)
        ex.set_ylim(1e-9, 1.0)
        ex.set_ylabel("|η − η_RCWA|")
        ex.set_xlabel(xl)
        ex.set_xlim(x[0], x[-1])
        if ex.get_legend_handles_labels()[0]:
            ex.legend(frameon=False, fontsize=6.8, loc="upper right", ncol=1, handlelength=2.2, labelspacing=0.2)
        s = stats(m)
        if TE:
            txt = f"rms error vs RCWA:  Kogelnik {s['kogelnik'][0]:.1e}   F-PSM {s['fpsm'][0]:.1e}"
        else:
            txt = (f"rms error vs RCWA (TM):  Kogelnik {s['kogelnik'][0]:.1e}  (relative, floor 10⁻³: {s['kogelnik_rel'][0]:.1e})\n"
                   "F-PSM: no TM formulae (2026-10-16)")
        ex.text(0.012, 0.06, txt, transform=ex.transAxes, fontsize=8.5, color=INK,
                bbox=dict(fc="#fcfcfb", ec="none", alpha=0.85, pad=1.5))
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=9, ncol=len(l) if len(l) <= 4 else 3, loc="lower center",
               bbox_to_anchor=(0.5, -0.035 if len(l) <= 4 else -0.075))
    axs[0, 0].set_title(f"Angle sweep at λ = {c['lam_nm']:.0f} nm", fontsize=9.5, color=MUTED, loc="left")
    axs[0, 1].set_title(f"Wavelength sweep at θ = {c['theta_bragg_deg']:.2f}°", fontsize=9.5, color=MUTED, loc="left")
    if PAPER:
        os.makedirs(os.path.join(FIG, "paper"), exist_ok=True)
        fig.savefig(os.path.join(FIG, "paper", f"compare_{name}{SUF}.svg"), bbox_inches="tight")
        fig.savefig(os.path.join(FIG, "paper", f"compare_{name}{SUF}.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(FIG, "paper", f"compare_{name}{SUF}.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)
    else:
        fig.suptitle(f"{c['title']}\n" + rf"{POL}, index-matched, $n_0$ = {c['n0']}, $n_1$ = {c['n1']}, d = {c['d_um']:g} µm, "
                     rf"slant Ψ = {c['psi'] / deg:.1f}°" + "      (grey band: outside F-PSM's stated domain)",
                     x=0.06, ha="left", fontsize=11, color=INK)
        fig.savefig(os.path.join(FIG, f"compare_{name}{SUF}.png"), dpi=170, bbox_inches="tight")
        plt.close(fig)

    for tag, x, m, z in (("angle", th, A, None if zm is None else zm[0]),
                         ("wavelength", lam, W, None if zm is None else zm[1])):
        s_all, ok = stats(m), m["in_domain"] & valid(m)
        s_in = stats(m, ok)
        zdev = np.nan
        if z is not None and len(z) == len(x):
            zdev = float(np.max(np.abs(z[:, 1] - m["kogelnik"])))
        row = [name, tag, len(x), int(ok.sum()), m["rcwa_conv"],
               s_all["kogelnik"][0], s_all["kogelnik"][1], s_all["fpsm"][0], s_all["fpsm"][1],
               s_in["kogelnik"][0], s_in["fpsm"][0], s_in["fpsm"][1], zdev]
        jj = 0 if tag == "angle" else 1
        for pre, *_ in ANSYS_RCWA:                 # Python RCWA, Kogelnik, F-PSM scored against each Ansys RCWA
            tab = None if ansys[pre] is None else ansys[pre][jj]
            if tab is None or not len(tab):
                row += [np.nan] * 6
                continue
            tab = tab[np.argsort(tab[:, 0])]
            good = np.isfinite(np.interp(tab[:, 0], x, np.where(valid(m), 0.0, np.nan)))
            d = {k: (np.interp(tab[:, 0], x, np.nan_to_num(m[k])) - tab[:, 1]) for k in ("rcwa", "kogelnik", "fpsm", "rcwa_epslin")}
            row += [float(np.sqrt(np.mean(d["rcwa"] ** 2))), float(np.max(np.abs(d["rcwa"]))),
                    float(np.sqrt(np.mean(d["kogelnik"][good] ** 2))),
                    float(np.sqrt(np.mean(d["fpsm"][good] ** 2))) if TE else np.nan,
                    float(np.sqrt(np.mean(d["rcwa_epslin"] ** 2))), float(np.max(np.abs(d["rcwa_epslin"])))]
        if not TE:
            row = [POL] + row + [s_all["kogelnik_rel"][0], s_all["kogelnik_rel"][1],
                                 float(np.sqrt(np.mean(((m["rcwa"] - m["rcwa_te"])[valid(m)]) ** 2))),
                                 float(np.max(np.abs((m["rcwa"] - m["rcwa_te"])[valid(m)])))]
        summary.append(row)


def main():
    os.makedirs(RES, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    summary = []
    for name in (sys.argv[1:] or CASES):
        one_case(name, CASES[name], summary)
        print("done", name, POL)
    head = ["case", "sweep", "n_points", "n_in_fpsm_domain", "rcwa_L10_vs_L14_max",
            "kogelnik_rms", "kogelnik_max", "fpsm_rms", "fpsm_max",
            "kogelnik_rms_in_domain", "fpsm_rms_in_domain", "fpsm_max_in_domain",
            "zemax_native_vs_my_kogelnik_max"]
    for pre, *_ in ANSYS_RCWA:
        head += [f"python_rcwa_vs_{pre}_rms", f"python_rcwa_vs_{pre}_max", f"kogelnik_vs_{pre}_rms", f"fpsm_vs_{pre}_rms",
                 f"python_rcwa_epslin_vs_{pre}_rms", f"python_rcwa_epslin_vs_{pre}_max"]
    if not TE:
        head = ["pol"] + head + ["kogelnik_rel_rms", "kogelnik_rel_max", "rcwa_TM_vs_TE_rms", "rcwa_TM_vs_TE_max"]
    if PAPER:
        return
    with open(os.path.join(RES, f"summary{SUF}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(head)
        for r in summary:
            w.writerow([f"{v:.3e}" if isinstance(v, float) else v for v in r])
    for r in summary:
        print(" | ".join(f"{v:.2e}" if isinstance(v, float) else str(v) for v in r))


if __name__ == "__main__":
    main()
