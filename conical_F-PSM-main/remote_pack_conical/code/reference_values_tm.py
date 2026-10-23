"""Reference values for the TM remote briefs: the converged Python referee (L = 14), TE and TM,
exact and eps-linear permittivity, plus Kogelnik TE/TM, at Bragg and at the discriminating
off-Bragg angles of every benchmark grating; and one absorbing point.  Writes
../evidence/tm_reference_values.md (markdown tables, quoted verbatim by the briefs).  2026-10-16."""
import os
import sys

import numpy as np

import fpsm
import kogelnik

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # Greek letters on a cp1252 console
import rcwa
from cases import CASES, deg

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OFF = {"lumerical_R": (0.0, -28.0, -14.0, 7.5), "lumerical_T": (0.0, 5.0), "zemax_kba_R": (-45.0, 0.0, -22.5),
       "unslanted_T": (33.8627180249, 36.0), "unslanted_T_rot5": (38.8627180249, 41.0)}
L = 14


def point(c, th, pol, form="exact", chi0=0.0, chi1=0.0):
    k0 = 2 * np.pi / c["lam_nm"]
    psi = c["psi"] if np.cos(th * deg + c["psi"]) >= 0 else c["psi"] + np.pi
    a = (k0, c["n0"], c["n1"], c["K"], psi, th * deg + psi, c["d_nm"])
    r = rcwa.solve(*a, L=L, eps_form=form, pol=pol, chi0=chi0, chi1=chi1)
    refl = int(fpsm.fpsm_detail(*a)["case"]) == -1                       # the note's own T/R rule (eq. 8), on the geometry
    first = r["DE_R"][L - 1] if refl else r["DE_T"][L - 1]
    t0 = r["DE_T"][L]
    kog = float(kogelnik.efficiency(*a, pol=pol, chi0=chi0, chi1=chi1))
    return first, t0, r["sum_DE"], kog, refl


lines = ["# TM reference values for the remote briefs\n",
         "*Computed 2026-10-16 by `code/reference_values_tm.py` with the gated referee `rcwa.py` (L = 14) and "
         "`kogelnik.py`; TE columns reproduce the 30 September benchmark's gate tables. Efficiencies per 1 W incident; "
         "angles inside the medium; R+1 / T+1 is the first order on the entry / far face. "
         "\"eps-linear\" is the OpticStudio DLL's permittivity expansion.*\n",
         "\n## Converged referee, both polarisations\n",
         "| case | θ (°) | λ (nm) | face | TE referee | TM referee | TM, eps-linear | TM/TE | TM T0 | Kogelnik TE | Kogelnik TM (κ cos 2θc) |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
for name, c in CASES.items():
    for th in OFF[name]:
        te, _, _, kte, refl = point(c, th, "TE")
        tm, t0, s, ktm, _ = point(c, th, "TM")
        tml, _, _, _, _ = point(c, th, "TM", "linear")
        assert abs(s - 1.0) < 1e-10, (name, th, s)
        lines.append(f"| `{name}` | {th:.4f} | {c['lam_nm']:.0f} | {'R+1' if refl else 'T+1'} | {te:.6f} | **{tm:.6f}** | {tml:.6f} | "
                     f"{tm / te:.4f} | {t0:.4f} | {kte:.6f} | {ktm:.6f} |")

lines += ["\n## The DLL probe of 2026-09-18 (10 layers per period, eps-linear) against the converged TM referee\n",
          "| θ (°) | DLL TE (10 layers) | DLL TM (10 layers) | DLL TM/TE | referee TM/TE (exact) | referee TM/TE (eps-linear) |",
          "|---|---|---|---|---|---|"]
dll = {0.0: (0.7965196, 0.7272546), -14.0: (0.0245959, 0.0245981), 7.5: (0.0978717, 0.052048)}
c = CASES["lumerical_R"]
for th, (dte, dtm) in dll.items():
    te, *_ = point(c, th, "TE")
    tm, *_ = point(c, th, "TM")
    tel, *_ = point(c, th, "TE", "linear")
    tml, *_ = point(c, th, "TM", "linear")
    lines.append(f"| {th:.1f} | {dte:.7f} | {dtm:.7f} | {dtm / dte:.4f} | {tm / te:.4f} | {tml / tel:.4f} |")

lines += ["\n## One absorbing point (for the optional Lumerical loss probe)\n",
          "G1 at Bragg (θ = 0°, 550 nm) with χ₀ = 1×10⁻³ (imaginary part of the mean index, n = 1.5 − 0.001 i), χ₁ = 0:\n",
          "| pol | first order R+1 | T0 | Σ over all orders (1 − absorbed fraction) |", "|---|---|---|---|"]
for pol in ("TE", "TM"):
    f, t0, s, _, _ = point(c, 0.0, pol, chi0=1e-3)
    lines.append(f"| {pol} | {f:.6f} | {t0:.6f} | {s:.6f} |")
lines.append("\nA solver whose imaginary index is silently ignored returns Σ = 1.000000 and the lossless first-order value; "
             "that is the gate of the optional step.\n")

out = os.path.join(ROOT, "evidence", "tm_reference_values.md")
with open(out, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print("->", out)
