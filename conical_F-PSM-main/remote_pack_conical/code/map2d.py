"""The Ansys Lumerical article's third figure -- first-order efficiency over
(replay angle, wavelength) -- for RCWA, Kogelnik and F-PSM, with the two error
maps against RCWA.   Usage: python map2d.py [case]   (default lumerical_R)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import fpsm
import kogelnik
import rcwa_te
from cases import CASES, deg

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
name = sys.argv[1] if len(sys.argv) > 1 else "lumerical_R"
c = CASES[name]
th = np.arange(c["theta_deg"][0], c["theta_deg"][-1] + 1e-9, 1.0)
lam = np.linspace(c["lam_sweep_nm"][0], c["lam_sweep_nm"][-1], 81)
TH, LAM = np.meshgrid(th, lam)
a = (2 * np.pi / LAM, c["n0"], c["n1"], c["K"], c["psi"], TH * deg + c["psi"], c["d_nm"])
det = fpsm.fpsm_detail(*a)
eR, eT = rcwa_te.signal_efficiency(*a, L=10)
refl = np.nanmax(eR) > np.nanmax(eT)
ref = eR if refl else eT
kog = kogelnik.efficiency(*a)
np.savez_compressed(os.path.join(ROOT, "results", f"map2d_{name}.npz"), theta_deg=th, lambda_nm=lam,
                    rcwa=ref, kogelnik=kog, fpsm=det["eta"], in_domain=det["in_domain"])

plt.rcParams.update({"font.size": 9, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb"})
fig, axs = plt.subplots(2, 3, figsize=(14, 7.8), sharex=True, sharey=True, gridspec_kw=dict(hspace=0.28))
for ax, (z, t) in zip(axs[0], ((ref, "RCWA (rigorous, L = 10)"), (kog, "Kogelnik"), (det["eta"], "F-PSM"))):
    im = ax.pcolormesh(th, lam, z, vmin=0, vmax=1, cmap="Blues", shading="nearest")
    ax.set_title(t, loc="left", fontsize=10)
fig.colorbar(im, ax=axs[0], pad=0.012, label=("R" if refl else "T") + "s power, order +1")
norm = LogNorm(1e-7, 1e-1)
axs[1, 0].axis("off")
NOTE = ("Top row: the three models on one\ncolour scale.\n\n"
        "Bottom row: absolute error against\nRCWA, log scale (darker = worse).\n\n"
        "Orange contour: edge of the working\ndomain F-PSM's note states for itself\n(R1-R6 / T1-T6).  The narrow strips\n"
        "are neighbour orders crossing cut-off\n(Rayleigh anomalies), where the 1/Γ\nsurface terms diverge.")
axs[1, 0].text(-0.12, 0.95, NOTE, transform=axs[1, 0].transAxes, va="top", fontsize=9, color="#52514e")
for ax, (z, t) in zip(axs[1, 1:], ((kog, "|Kogelnik − RCWA|"), (det["eta"], "|F-PSM − RCWA|"))):
    e = np.maximum(np.abs(z - ref), 1e-7)
    im2 = ax.pcolormesh(th, lam, e, norm=norm, cmap="Blues", shading="nearest")
    rms = float(np.sqrt(np.nanmean((z - ref) ** 2)))
    ax.set_title(f"{t}\nrms {rms:.1e}   max {np.nanmax(np.abs(z - ref)):.1e}", loc="left", fontsize=10)
    ax.set_xlabel("Replay angle in the medium, θ [°]")
axs[1, 2].contour(th, lam, det["in_domain"].astype(float), levels=[0.5], colors="#eb6834", linewidths=1.2)
axs[1, 2].plot([], [], color="#eb6834", lw=1.2, label="boundary of F-PSM's stated domain")
axs[1, 2].legend(frameon=False, loc="lower left", fontsize=8.5)
fig.colorbar(im2, ax=axs[1], pad=0.012, label="absolute efficiency error")
axs[0, 0].set_ylabel("Replay wavelength [nm]")
axs[1, 1].set_ylabel("Replay wavelength [nm]")
axs[1, 1].tick_params(labelleft=True)
fig.suptitle(f"{c['title']} — TE, $n_0$ = {c['n0']}, $n_1$ = {c['n1']}, d = {c['d_um']:g} µm", x=0.06, ha="left", fontsize=11)
fig.savefig(os.path.join(ROOT, "figures", f"map2d_{name}.png"), dpi=160, bbox_inches="tight")
print("rms kogelnik", float(np.sqrt(np.nanmean((kog - ref) ** 2))), "rms fpsm", float(np.sqrt(np.nanmean((det['eta'] - ref) ** 2))),
      "fpsm in-domain rms", float(np.sqrt(np.nanmean(((det['eta'] - ref)[det['in_domain']]) ** 2))), "in-domain fraction", float(det["in_domain"].mean()))
