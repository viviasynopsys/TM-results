"""figures_conical -- the two figures of the report from results/sweep_504.csv and
results/reference_conical.csv.

    python figures_conical.py
    -> figures/sweep_errors_vs_mu.png     (the report's Fig. 2 regenerated: max error vs mu per family)
    -> figures/reference_conical_bench.png (referee s and p curves on the conical Bragg cones of the
                                            four benchmark gratings, corrected model dotted)
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES, FIG = os.path.join(ROOT, "results"), os.path.join(ROOT, "figures")


def rd(name):
    with open(os.path.join(RES, name), newline="") as f:
        return list(csv.DictReader(f))


def fig_sweep():
    rows = rd("sweep_504.csv")
    fams = ["Slanted T", "Slanted R", "Unslanted T"]
    mus = sorted({float(r["mu"]) for r in rows})
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for ax, fam in zip(axes, fams):
        sub = [r for r in rows if r["family"] == fam]
        L = [max(float(r["err_leading"]) for r in sub if abs(float(r["mu"]) - mu) < 1e-9) for mu in mus]
        C = [max(float(r["err_corrected"]) for r in sub if abs(float(r["mu"]) - mu) < 1e-9) for mu in mus]
        Cc = [max([float(r["err_corrected"]) for r in sub if abs(float(r["mu"]) - mu) < 1e-9 and r["neighbours_closed"] == "True"] or [np.nan]) for mu in mus]
        ax.loglog(mus, L, "o--", label="leading")
        ax.loglog(mus, C, "o-", label="cubic + faces")
        ax.loglog(mus, Cc, "s:", color="k", label="closed neighbours")
        ax.set_title(fam)
        ax.set_xlabel("mu = n1/n0")
        ax.set_xticks(mus)
        ax.set_xticklabels([f"{m:.3f}" for m in mus], fontsize=7)
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("max all-polarisation efficiency error")
    axes[0].legend(fontsize=8)
    fig.suptitle("504-point sweep of the report, re-run here: closed forms vs the independent conical referee (19 orders)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "sweep_errors_vs_mu.png"), dpi=140)
    plt.close(fig)


def fig_reference():
    rows = [r for r in rd("reference_conical.csv") if r["set"] == "BENCH"]
    cases = ["lumerical_R", "lumerical_T", "zemax_kba_R", "unslanted_T"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, name in zip(axes.ravel(), cases):
        sub = [r for r in rows if r["case"] == name]
        phis = sorted({float(r["phi_deg"]) for r in sub})
        for i, phi in enumerate(phis):
            s = [r for r in sub if float(r["phi_deg"]) == phi]
            th = np.array([float(r["theta_deg"]) for r in s])
            es, ep = np.array([float(r["eta_s"]) for r in s]), np.array([float(r["eta_p"]) for r in s])
            cs = np.array([float(r["eta_s_corr"]) if r["eta_s_corr"] not in ("", None) else np.nan for r in s])
            cp = np.array([float(r["eta_p_corr"]) if r["eta_p_corr"] not in ("", None) else np.nan for r in s])
            ax.plot(th, es, "-", color=f"C{i}", label=f"phi {phi:.0f}: s")
            ax.plot(th, ep, "--", color=f"C{i}", label=f"phi {phi:.0f}: p")
            ax.plot(th, cs, ":", color="k", lw=0.8)
            ax.plot(th, cp, ":", color="k", lw=0.8)
        ax.set_title(name)
        ax.set_xlabel("theta (deg)")
        ax.set_ylabel("order-1 efficiency")
        ax.legend(fontsize=7, ncol=2)
        ax.grid(alpha=0.3)
    fig.suptitle("Conical Bragg cones of the benchmark gratings: referee (colour), corrected closed form (black dotted)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "reference_conical_bench.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_sweep()
    fig_reference()
    print("-> figures/sweep_errors_vs_mu.png, figures/reference_conical_bench.png")
