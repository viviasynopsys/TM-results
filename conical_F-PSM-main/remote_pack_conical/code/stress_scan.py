"""Where does F-PSM fail?  Test over a wider range than the note's own
135 + 135 gratings -- high modulation, thickness, slant.  TE or TM (2026-10-16).

    python stress_scan.py                 TE: the 30 September scan (F-PSM, Kogelnik, referee)
    python stress_scan.py --pol TM        TM: referee and Kogelnik only -- F-PSM has no TM formulae
                                          (fpsm.TM_ABSENT); this is the reference set a TM sheet
                                          will be scored against
    ... --figures-only                    redraw from the CSV

Construction follows the note's Section 6 realistic-units cases (and the paper's
Section 8 families): n0 = 1.5, recording wavelength 500 nm, one fixed reference
beam at internal incidence 33.8627180249 deg, the slant swept, the grating
replayed by its own reference beam (theta_r = theta_c = theta_inc + psi,
K = 4 pi n0 cos(theta_r) / lam_r, signed).  Metric = the note's: rms efficiency error
against RCWA across the 400-600 nm band (101 points), lossless.

    reflection family   psi = -60 .. +30 deg (step 5, psi = 0 excluded: the note
                        says the unslanted reflection case is a separate problem)
    transmission family psi =  50 .. 130 deg (step 5; 90 = unslanted)
    n1 in {0.005, 0.01, 0.02, 0.03, 0.045, 0.06, 0.08, 0.12, 0.20}
    d  in {2, 5, 10, 20, 50} um

Output: ../results/stress_scan<SUF>.csv, ../figures/stress_*<SUF>.png  (SUF = "" TE, "_TM" TM).
The TE path is the 30 September script unchanged (geometry, ladder, columns, figures); TM rows add
``pol`` and ``kogelnik_rel_rms`` (|dK| / max(eta_ref, 1e-3)) and leave every F-PSM column empty.
"""
import csv
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import fpsm
import kogelnik
import rcwa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
deg = np.pi / 180
N0, LAM_R, TH_INC = 1.5, 500.0, 33.8627180249
LAMS = np.linspace(400.0, 600.0, 101)
N1S = [0.005, 0.01, 0.02, 0.03, 0.045, 0.06, 0.08, 0.12, 0.20]
DS_UM = [2, 5, 10, 20, 50]
FAMILIES = {"reflection": [p for p in range(-60, 31, 5) if p != 0],
            "transmission": list(range(50, 131, 5))}
L_LADDER = (10, 16, 24, 36, 54)   # RCWA truncation; see converged_L
CONV_TOL = 1.0e-8
REL_FLOOR = 1e-3


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


def converged_L(a, k0_probe):
    """Smallest truncation on the ladder whose signal efficiencies (both faces)
    agree with the next rung to CONV_TOL on an 11-point probe of the band.
    L = 10 is NOT enough for the long-period, many-order gratings near
    psi ~ 50-60 deg (3.8e-2 off at psi = 55, n1 = 0.2: VERIFIED 2026-09-17), so
    the referee must earn its own trust point by point.  Returns (L, residual);
    L = -1 if the ladder is exhausted -- the caller then records the grating as
    unrefereed instead of quoting an error against an unconverged reference."""
    probe = (k0_probe,) + a[1:]
    prev = rcwa.signal_efficiency(*probe, L=L_LADDER[0], pol=POL)
    resid = np.inf
    for L_lo, L_hi in zip(L_LADDER[:-1], L_LADDER[1:]):
        nxt = rcwa.signal_efficiency(*probe, L=L_hi, pol=POL)
        resid = float(max(np.max(np.abs(prev[0] - nxt[0])), np.max(np.abs(prev[1] - nxt[1]))))
        if resid < CONV_TOL:
            return L_lo, resid
        prev = nxt
    return -1, resid


def run():
    rows = []
    t0 = time.time()
    k0 = 2 * np.pi / LAMS
    for fam, psis in FAMILIES.items():
        for psi_d in psis:
            psi = psi_d * deg
            thc = TH_INC * deg + psi
            K = fpsm.K_from_recording(N0, LAM_R, thc)
            for n1 in N1S:
                for d_um in DS_UM:
                    a = (k0, N0, n1, K, psi, thc, d_um * 1e3)
                    det = fpsm.fpsm_detail(*a)                    # TE lossless: geometry, flags, TE F-PSM
                    ok = np.isfinite(det["eta"])
                    L_used, resid = converged_L(a, k0[::10])
                    if L_used < 0:
                        row = dict(family=fam, psi_deg=psi_d, n1=n1, d_um=d_um, geometry_at_bragg="unrefereed",
                                   rcwa_L=-1, rcwa_residual=resid)
                        if not TE:
                            row["pol"] = POL
                        rows.append(row)
                        continue
                    eR, eT = rcwa.signal_efficiency(*a, L=L_used, pol=POL)
                    # the face is decided at Bragg (lam_r) by the note's own rule, eq. (8)
                    bragg = fpsm.fpsm_detail(2 * np.pi / LAM_R, N0, n1, K, psi, thc, d_um * 1e3)
                    case_b = int(bragg["case"])
                    ref = eR if case_b == -1 else eT
                    kog = kogelnik.efficiency(*a, pol=POL)
                    same = ok & (det["case"] == case_b)
                    fe = det["eta"] if TE else np.full_like(ref, np.nan)
                    ef = (fe - ref)[same]
                    ek = (kog - ref)[same]
                    ind = same & det["in_domain"]
                    efi, eki = (fe - ref)[ind], (kog - ref)[ind]
                    row = dict(
                        family=fam, psi_deg=psi_d, n1=n1, d_um=d_um,
                        geometry_at_bragg={-1: "reflection", 1: "transmission", 0: "closed"}[case_b],
                        rcwa_L=L_used, rcwa_residual=resid,
                        n_points=int(same.sum()), n_closed=int((~ok).sum()),
                        n_flag_flips=int((ok & (det["case"] != case_b)).sum()),
                        n_in_domain=int((det["in_domain"] & same).sum()),
                        bragg_in_domain=int(bool(bragg["in_domain"])),
                        bragg_min_gamma_over_beta=float(min(bragg["gamma0"], bragg["gamma1"]) / (2 * np.pi / LAM_R * N0)),
                        fpsm_rms_in_domain_pts=float(np.sqrt(np.mean(efi ** 2))) if efi.size and TE else np.nan,
                        fpsm_max_in_domain_pts=float(np.max(np.abs(efi))) if efi.size and TE else np.nan,
                        kogelnik_rms_in_domain_pts=float(np.sqrt(np.mean(eki ** 2))) if eki.size else np.nan,
                        fpsm_rms=float(np.sqrt(np.mean(ef ** 2))) if ef.size and TE else np.nan,
                        fpsm_max=float(np.max(np.abs(ef))) if ef.size and TE else np.nan,
                        kogelnik_rms=float(np.sqrt(np.mean(ek ** 2))) if ek.size else np.nan,
                        fpsm_out_of_unit=float(np.nanmax(np.r_[0.0, (fe[same] - 1.0), -fe[same]])) if TE else np.nan,
                        domain_fail_or=int(np.bitwise_or.reduce(det["domain_fail"][same])) if same.any() else -1,
                        peak_eta_rcwa=float(ref[same].max()) if same.any() else np.nan)
                    if not TE:
                        rel = (np.abs(kog - ref) / np.maximum(ref, REL_FLOOR))[same]
                        row["pol"] = POL
                        row["kogelnik_rel_rms"] = float(np.sqrt(np.mean(rel ** 2))) if rel.size else np.nan
                        row["kogelnik_rel_max"] = float(np.max(rel)) if rel.size else np.nan
                    rows.append(row)
        print(f"{fam} ({POL}): done in {time.time() - t0:.0f} s", flush=True)
    return rows


def save(rows):
    p = os.path.join(ROOT, "results", f"stress_scan{SUF}.csv")
    with open(p, "w", newline="") as f:
        full = max(rows, key=len)
        w = csv.DictWriter(f, fieldnames=list(full), restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4e}" if isinstance(v, float) else v) for k, v in r.items()})
    return p


def figures(rows):
    rows = [r for r in rows if r["geometry_at_bragg"] != "unrefereed"]
    plt.rcParams.update({"font.size": 8.5, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                         "savefig.facecolor": "#fcfcfb", "axes.spines.top": False, "axes.spines.right": False})
    norm = LogNorm(1e-7, 1e-1)
    series = (("fpsm_rms", "F-PSM"), ("kogelnik_rms", "Kogelnik")) if TE else \
             (("kogelnik_rms", "Kogelnik (TM)"), ("peak_eta_rcwa", "referee peak η (TM)"))
    for fam, psis in FAMILIES.items():
        fig, axs = plt.subplots(2, len(DS_UM), figsize=(15, 5.6), sharex=True, sharey=True)
        for r_, (key, lab) in enumerate(series):
            for c_, d_um in enumerate(DS_UM):
                G = np.full((len(N1S), len(psis)), np.nan)
                D = np.zeros_like(G, dtype=bool)
                for r in rows:
                    if r["family"] == fam and r["d_um"] == d_um:
                        i, j = N1S.index(r["n1"]), psis.index(r["psi_deg"])
                        G[i, j] = max(r[key], 1e-7) if np.isfinite(r[key]) else np.nan
                        D[i, j] = bool(r["bragg_in_domain"])
                ax = axs[r_, c_]
                if key == "peak_eta_rcwa":
                    im2 = ax.pcolormesh(np.arange(len(psis) + 1), np.arange(len(N1S) + 1), G, vmin=0, vmax=1,
                                        cmap="Blues", edgecolors="#fcfcfb", linewidth=1.0)
                else:
                    im = ax.pcolormesh(np.arange(len(psis) + 1), np.arange(len(N1S) + 1), G, norm=norm,
                                       cmap="Blues", edgecolors="#fcfcfb", linewidth=1.0)
                if key == "fpsm_rms" or (not TE and r_ == 0):
                    jj, ii = np.meshgrid(np.arange(len(psis)), np.arange(len(N1S)))
                    ax.plot(jj[D] + 0.5, ii[D] + 0.5, ls="none", marker="o", ms=2.6, color="#eb6834")
                ax.set_xticks(np.arange(len(psis))[::3] + 0.5, [str(p) for p in psis[::3]])
                ax.set_yticks(np.arange(len(N1S)) + 0.5, [f"{v:g}" for v in N1S])
                if r_ == 0:
                    ax.set_title(f"d = {d_um} µm", fontsize=9.5, loc="left", color="#52514e")
                if r_ == 1:
                    ax.set_xlabel("slant Ψ [°]")
                if c_ == 0:
                    ax.set_ylabel(f"{lab}\nmodulation $n_1$")
        if TE:
            cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.015)
            cb.set_label("rms efficiency error vs RCWA, 400–600 nm")
            fig.suptitle(f"{fam.capitalize()} family — rms error against RCWA (TE, lossless, $n_0$ = 1.5, reference beam 33.86° internal).  "
                         "Orange dot: the Bragg point (500 nm) is inside F-PSM's stated domain.", x=0.06, ha="left", fontsize=10.5)
        else:
            cb = fig.colorbar(im, ax=axs[0], shrink=0.9, pad=0.015)
            cb.set_label("Kogelnik rms error vs RCWA (TM)")
            cb2 = fig.colorbar(im2, ax=axs[1], shrink=0.9, pad=0.015)
            cb2.set_label("referee peak η over 400–600 nm (TM)")
            fig.suptitle(f"{fam.capitalize()} family — TM, lossless, $n_0$ = 1.5, reference beam 33.86° internal.  "
                         "Top: Kogelnik (κ cos 2θ_c) rms error against the TM referee; bottom: the TM referee's peak efficiency.  "
                         "Orange dot: Bragg point inside F-PSM's TE domain.  F-PSM: no TM formulae (2026-10-16).",
                         x=0.06, ha="left", fontsize=10)
        fig.savefig(os.path.join(ROOT, "figures", f"stress_{fam}{SUF}.png"), dpi=160, bbox_inches="tight")
        plt.close(fig)

    if not TE:
        return
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.2), sharex=True, sharey=True)
    for ax, fam in zip(axs, FAMILIES):
        rr = [r for r in rows if r["family"] == fam and np.isfinite(r["fpsm_rms"])]
        ind = np.array([bool(r["bragg_in_domain"]) for r in rr])
        x = np.array([max(r["kogelnik_rms"], 1e-9) for r in rr])
        y = np.array([max(r["fpsm_rms"], 1e-9) for r in rr])
        ax.loglog(x[~ind], y[~ind], ls="none", marker="o", ms=4, mfc="none", mec="#2a78d6", mew=1.0,
                  label=f"Bragg point outside stated domain ({(~ind).sum()})")
        ax.loglog(x[ind], y[ind], ls="none", marker="o", ms=4.5, color="#eb6834", label=f"Bragg point inside stated domain ({ind.sum()})")
        ax.plot([1e-8, 1], [1e-8, 1], color="#52514e", lw=1)
        ax.text(3e-6, 5e-6, "F-PSM = Kogelnik", rotation=40, color="#52514e", fontsize=8, ha="left", va="bottom")
        ax.set_xlim(1e-6, 1)
        ax.set_ylim(1e-8, 1)
        ax.grid(True, color="#e4e3df", lw=0.6)
        ax.set_title(f"{fam.capitalize()} family   ({int(np.sum(y > x))} of {len(rr)} gratings above the line: F-PSM worse)",
                     fontsize=9.5, loc="left", color="#52514e")
        ax.set_xlabel("Kogelnik rms error vs RCWA")
        ax.legend(frameon=False, loc="lower right", fontsize=8.5)
    axs[0].set_ylabel("F-PSM rms error vs RCWA")
    fig.suptitle("Is F-PSM ever worse than Kogelnik?  One point per grating; rms over 400–600 nm.", x=0.06, ha="left", fontsize=10.5)
    fig.savefig(os.path.join(ROOT, "figures", "stress_fpsm_vs_kogelnik.png"), dpi=160, bbox_inches="tight")
    plt.close(fig)


def load():
    """Re-read ../results/stress_scan<SUF>.csv so figures can be redrawn without the ~15 min scan."""
    ints = {"psi_deg", "d_um", "n_points", "n_closed", "n_flag_flips", "n_in_domain", "bragg_in_domain",
            "domain_fail_or", "rcwa_L"}
    rows = []
    for r in csv.DictReader(open(os.path.join(ROOT, "results", f"stress_scan{SUF}.csv"))):
        rows.append({k: (v if k in ("family", "geometry_at_bragg", "pol") else np.nan if v == "" else
                         int(float(v)) if k in ints else float(v)) for k, v in r.items()})
    return rows


if __name__ == "__main__":
    if "--figures-only" in sys.argv:
        figures(load())
    else:
        rows = run()
        print(save(rows))
        figures(rows)
