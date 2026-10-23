"""One picture for the whole benchmark: rms error of every model against the rigorous referee, per
grating and sweep, on a log axis.  Reads ../results/summary<SUF>.csv (written by compare.py) and
writes ../figures/summary_rms<SUF>.png.   ``--pol TM`` selects the TM tables (2026-10-16); in TM
the closed-form facet holds Kogelnik only, because F-PSM has no TM formulae (fpsm.TM_ABSENT).

Two facets sharing the y rows and the x axis: the closed-form theories (Kogelnik, F-PSM) and the
two Ansys RCWA solvers (OpticStudio DLL, Lumerical).  Colours are the slots compare.py uses for
the same entities (blue, orange; yellow, magenta -- the Ansys pair drawn with the palette's darker
steps as marker edges).  Each facet's pair was validated with the data-viz palette checker
(all-pairs mode, light surface); identity is carried by colour AND marker shape, and a legend is
always present.  The hollow yellow square is the DLL scored against the referee with the DLL's own
permittivity expansion (eps = n0^2 + 2 n0 dn cos), the diagnostic that explains its distance from
the exact referee.  The full table behind the picture is results/summary<SUF>.csv.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
C_KOG, C_FPSM = "#2a78d6", "#eb6834"
C_ZRC, C_ZRC_EDGE, C_LUM, C_LUM_EDGE = "#eda100", "#c98500", "#e87ba4", "#d55181"

TITLES = {"lumerical_R": "Reflection VHG, Lumerical article (n₁ 0.04, 6 µm, Ψ 14°)",
          "zemax_kba_R": "Reflection VHG, Zemax article (n₁ 0.05, 5 µm, Ψ 22.5°)",
          "lumerical_T": "Transmission variant, Ψ 76° (outside F-PSM's domain)",
          "unslanted_T": "Unslanted transmission, the note's T-a (n₁ 0.03, 10 µm)",
          "unslanted_T_rot5": "T-a rotated 5°, the DLL substitute (Ψ −95°)"}
ORDER = ["lumerical_R", "zemax_kba_R", "lumerical_T", "unslanted_T", "unslanted_T_rot5"]
FLOOR = 1e-13
GATE = 1e-3


def _pol(argv):
    if "--pol" in argv:
        i = argv.index("--pol")
        pol = argv[i + 1].upper()
        del argv[i:i + 2]
        return pol
    return "TE"


POL = _pol(sys.argv)
TE = POL == "TE"
SUF = "" if TE else "_TM"


def val(r, col):
    try:
        v = float(r.get(col, "nan"))
    except ValueError:
        return np.nan
    return v


def main():
    paper = os.environ.get("COMPARE_PAPER", "") not in ("", "0")
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", f"summary{SUF}.csv"))))
    by = {(r["case"], r["sweep"]): r for r in rows}
    rows_y = []                                            # (y, case, sweep)
    y = 0.0
    for case in ORDER:
        for sweep in ("angle", "wavelength"):
            if (case, sweep) in by:
                rows_y.append((y, case, sweep))
                y += 1.0
        y += 0.7
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                         "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
                         "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF})
    if paper:
        plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42})
    fig, axs = plt.subplots(1, 2, figsize=(7.06, 3.7) if paper else (12.5, 6.2), sharey=True,
                            gridspec_kw=dict(wspace=0.06))
    closed = [("kogelnik_rms", "Kogelnik" if TE else "Kogelnik, TM (κ cos 2θ_c)", C_KOG, C_KOG, "o", 7.5)]
    if TE:
        closed.append(("fpsm_rms", "F-PSM", C_FPSM, C_FPSM, "D", 6.5))
    facets = [
        ("Closed-form theories" if TE else "Closed-form theory (F-PSM: no TM formulae, 2026-10-16)", closed),
        ("Ansys RCWA solvers", [("python_rcwa_vs_zemax_rcwa_rms", "OpticStudio RCWA DLL", C_ZRC, C_ZRC_EDGE, "s", 7.5),
                                ("python_rcwa_epslin_vs_zemax_rcwa_rms", "DLL vs referee, ε linear in n₁", SURF, C_ZRC_EDGE, "s", 7.5),
                                ("python_rcwa_vs_lumerical_rcwa_rms", "Lumerical RCWA", C_LUM, C_LUM_EDGE, "^", 8.5)]),
    ]
    for ax, (title, series) in zip(axs, facets):
        ax.grid(True, axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for yy, _, _ in rows_y:
            ax.axhline(yy, color=GRID, lw=0.6, zorder=0)
        ax.axvline(GATE, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
        drawn = 0
        for col, lab, fc, ec, mk, ms in series:
            xs, yy = [], []
            for y0, case, sweep in rows_y:
                v = val(by[(case, sweep)], col)
                if np.isfinite(v):
                    xs.append(max(v, FLOOR))
                    yy.append(y0)
            drawn += len(xs)
            ax.plot(xs, yy, ls="none", marker=mk, ms=ms, mfc=fc, mec=ec, mew=1.7, label=lab, zorder=3)
        if drawn == 0:
            ax.text(1e-7, rows_y[len(rows_y) // 2][0], f"no {POL} Ansys tables yet", color=MUTED, fontsize=9, ha="center", va="center")
        ax.set_xscale("log")
        ax.set_xlim(3e-14, 1e-1)
        ax.set_xticks([1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2])
        ax.set_title(title, loc="left", fontsize=10, color=MUTED)
        ax.tick_params(axis="y", length=0)
        ax.legend(frameon=False, fontsize=8 if paper else 8, loc="upper left", bbox_to_anchor=(0.0, -0.12 if paper else -0.09),
                  handletextpad=0.5, labelspacing=0.35)
        ax.text(GATE, rows_y[0][0] - 0.75, "10⁻³ threshold" if paper else "gate 10⁻³", color=MUTED, fontsize=7.5 if not paper else 8, ha="center", va="bottom")
        ax.set_xlabel("rms |η − η_referee| over the sweep", color=INK, fontsize=8.5 if not paper else 8)
    axs[0].set_yticks([yy for yy, _, _ in rows_y])
    short = {"lumerical_R": "G1", "zemax_kba_R": "G2", "lumerical_T": "G3", "unslanted_T": "G4", "unslanted_T_rot5": "G5"}
    if paper:
        axs[0].set_yticklabels([f"{short[c]} {s}" for _, c, s in rows_y], fontsize=8)
        if TE and len(rows_y) > 6:
            axs[1].text(3e-11, rows_y[6][0], "Lumerical 10⁻¹²;  DLL: no result (K_z = 0)", fontsize=7.5, color=MUTED, va="center", ha="left")
    else:
        axs[0].set_yticklabels([f"{TITLES[c]}\n{s} sweep" for _, c, s in rows_y], fontsize=7.8)
        if TE and len(rows_y) > 6:
            axs[1].text(3e-11, rows_y[6][0], "◄ Lumerical = referee to 10⁻¹²;  DLL: no result (K_z = 0 crashes the server)",
                        fontsize=7.5, color=MUTED, va="center", ha="left")
    axs[0].invert_yaxis()
    if not paper:
        fig.text(0.5, -0.06, f"Referee: Python RCWA ({POL}) with ε = (n₀ + n₁cos)², L = 10 (converged against L = 14).  "
                 f"Every value is in results/summary{SUF}.csv; per-point curves in results/compare_<case>{SUF}_*.csv.",
                 fontsize=8, color=MUTED, ha="center")
    if paper:
        os.makedirs(os.path.join(ROOT, "figures", "paper"), exist_ok=True)
        for ext in ("svg", "png", "pdf"):
            out = os.path.join(ROOT, "figures", "paper", f"summary_rms{SUF}.{ext}")
            fig.savefig(out, dpi=200, bbox_inches="tight")
            print("->", out)
        return
    if TE:
        fig.suptitle("Every model against the rigorous referee — five gratings, both sweeps\n"
                     "Ansys Zemax OpticStudio 2024 R1 RCWA DLL and Ansys Lumerical 2026 R1 RCWA; VERIFIED 2026-09-18, rescored 2026-09-30",
                     x=0.01, ha="left", fontsize=10.5, color=INK)
    else:
        fig.suptitle("TM: every model against the rigorous TM referee — five gratings, both sweeps\n"
                     "Ansys TM tables appear here when the remote runs return (remote_pack/).  Python legs VERIFIED 2026-10-16",
                     x=0.01, ha="left", fontsize=10.5, color=INK)
    out = os.path.join(ROOT, "figures", f"summary_rms{SUF}.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print("->", out)


if __name__ == "__main__":
    main()
