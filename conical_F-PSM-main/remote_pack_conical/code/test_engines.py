"""Referee checks: rcwa_te vs (1) the 'rigorous' values printed in the
F-PSM note, Section 6 -- an oracle independent of this code -- (2) energy
conservation, (3) an optional second RCWA/Kogelnik package (env PSMGRATING_PATH);
and kogelnik.py vs that package's Kogelnik.  Exit code 1 on any failure."""
import sys
import numpy as np
import fpsm, rcwa_te, kogelnik

deg = np.pi / 180
fails = []

def check(name, got, want, atol):
    ok = abs(got - want) <= atol
    print(f"{'PASS' if ok else 'FAIL'}  {name:40s} got {got:.10f}  want {want:.10f}  |diff| {abs(got-want):.1e}")
    if not ok: fails.append(name)

k5 = 2 * np.pi / 500
cases = [  # tag, args (k0,n0,n1,K,psi,thc,d), the note's rigorous eta, face
    ("T-a", (k5, 1.5, 0.03, -2.1006130130e-2, 90*deg, 123.8627180249*deg, 1e4), 0.5860347818, "T"),
    ("T-b", (k5, 1.5, 0.03, -1.5251047602e-2, 80*deg, 113.8627180249*deg, 5e3), 0.7518293859, "T"),
    ("R-a", (k5, 1.5, 0.03, 3.5674516079e-2, -15*deg, 18.8627180249*deg, 1e4), 0.9382471792, "R"),
    ("R-b", (k5, 1.5, 0.03, 3.7248996723e-2, -25*deg, 8.8627180249*deg, 5e3), 0.6142644669, "R"),
    ("R-dimless-1", (1/np.cos(25*deg), 1.0, 0.02, 2.0, 15*deg, 25*deg, 39.3593963), 0.2134538050, "R"),
    ("R-dimless-2", (1/np.cos(25*deg), 1.0, 0.02, 2.0, 15*deg, 25*deg, 118.0781889), 0.8190904574, "R"),
    ("R-dimless-3", (1/np.cos(25*deg), 1.0, 0.02, 2.0, 15*deg, 25*deg, 236.1563778), 0.9899977473, "R"),
]
for tag, a, want, face in cases:
    for L in (10,):
        r = rcwa_te.solve(*a, L=L)
        got = r["DE_R"][L-1] if face == "R" else r["DE_T"][L-1]
        check(f"RCWA L={L} vs note's rigorous {tag}", got, want, 2e-9)
        check(f"RCWA energy {tag}", r["sum_DE"], 1.0, 1e-10)
    case = int(fpsm.fpsm_detail(*a)["case"])
    check(f"FPSM case flag {tag}", case, 1 if face == "T" else -1, 0)

try:
    # Optional cross-check against a second, separately validated package: point
    # PSMGRATING_PATH at the folder that holds ``psmgrating/``.  Absent -> skipped.
    import os
    if not os.environ.get("PSMGRATING_PATH"):
        raise ImportError("PSMGRATING_PATH not set")
    sys.path.insert(0, os.environ["PSMGRATING_PATH"])
    from psmgrating import rcwa as xr, twowave as tw
    lam_r = 550e-9; n0, n1, d = 1.5, 0.04, 6e-6; psi = 14*deg
    K = fpsm.K_from_recording(n0, 550.0, 14*deg)
    for th in (-35, -28, -14, 0, 7.5):
        thc = (th + 14) * deg
        x = xr.solve(lam_c=550e-9, lam_r=lam_r, n0=n0, n1=n1, chi0=0, chi1=0, theta_c=thc, theta_r=14*deg, psi=psi, d=d, L=10)
        mine = rcwa_te.solve(2*np.pi/550, n0, n1, K, psi, thc, 6000.0, L=10)["DE_R"][9]
        check(f"RCWA vs explorer, Lumerical case th={th}", mine, x["DE_R"][9], 1e-10)
        kx = float(tw.efficiency("kogelnik", "sigma", n0=n0, n1=n1, d=d, lam_r=lam_r, theta_r=14*deg, lam_c=550e-9, theta_c=thc, psi=psi))
        check(f"Kogelnik vs explorer th={th}", float(kogelnik.efficiency(2*np.pi/550, n0, n1, K, psi, thc, 6000.0)), kx, 1e-9)
except ImportError as e:
    print("SKIP explorer cross-check (package not importable):", e)

print(f"\n{len(fails)} failure(s)"); sys.exit(1 if fails else 0)
