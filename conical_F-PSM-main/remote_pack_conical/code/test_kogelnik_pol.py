"""Gates for kogelnik.py (TE/TM, lossless/lossy).  Exit code 1 on any failure.  2026-10-16.

  K1  TE lossless == kogelnik_te.efficiency (the 30 September engine) bit for bit, every grating sweep
  K2  the general 2x2 path on TE lossless == K1 path to 1e-12 (validates the general solution)
  K3  general path == explorer psmgrating.twowave Kogelnik: TE lossy, TM lossless, TM lossy, to 1e-10
  K4  TM at theta_c = 45 deg (cos 2 theta_c = 0) gives zero coupling -> eta = 0
"""
import os
import sys

import numpy as np

import kogelnik
import kogelnik_te
from cases import CASES, deg

fails = []
skips = []


def check(name, got, want, atol):
    got, want = np.asarray(got, float), np.asarray(want, float)
    worst = float(np.max(np.abs(got - want)))
    ok = worst <= atol
    print(f"{'PASS' if ok else 'FAIL'}  {name:62s} max|diff| {worst:.2e}  (tol {atol:.0e})")
    if not ok:
        fails.append(name)


def sweep_args(c):
    th = c["theta_deg"]
    k0 = 2 * np.pi / c["lam_nm"]
    psi = np.where(np.cos(th * deg + c["psi"]) >= 0, c["psi"], c["psi"] + np.pi)
    return (k0, c["n0"], c["n1"], c["K"], psi, th * deg + psi, c["d_nm"])


print("K1  TE lossless against kogelnik_te (bit-identical)")
for name, c in CASES.items():
    a = sweep_args(c)
    old, new = kogelnik_te.efficiency(*a), kogelnik.efficiency(*a)
    ident = np.array_equal(old, new)
    print(f"{'PASS' if ident else 'FAIL'}  {'K1 ' + name:62s} bit-identical {ident}  max|diff| {np.max(np.abs(old - new)):.1e}")
    if not ident:
        fails.append("K1 " + name)

print("K2  general 2x2 path on TE lossless against the closed formulae")
for name, c in CASES.items():
    a = sweep_args(c)
    gen = kogelnik._general(*[np.asarray(v, float) for v in a], "TE", np.zeros_like(a[0]), np.zeros_like(a[0]))
    check(f"K2 {name}", gen, kogelnik_te.efficiency(*a), 1e-12)
    # wavelength sweep too
    lam = c["lam_sweep_nm"]
    k0 = 2 * np.pi / lam
    thb = c["theta_bragg_deg"]
    psi = c["psi"] if np.cos(thb * deg + c["psi"]) >= 0 else c["psi"] + np.pi
    a2 = (k0, c["n0"], c["n1"], c["K"], psi, thb * deg + psi, c["d_nm"])
    gen2 = kogelnik._general(*[np.asarray(v, float) for v in np.broadcast_arrays(*[np.asarray(v, float) for v in a2])],
                             "TE", np.zeros_like(k0), np.zeros_like(k0))
    check(f"K2 {name} wavelength", gen2, kogelnik_te.efficiency(*a2), 1e-12)

print("K3  against the explorer package's Kogelnik (TE lossy, TM lossless, TM lossy)")
xp = os.environ.get("PSMGRATING_PATH", r"C:\Users\vivia\Repos\psm-explorer\python")
try:
    sys.path.insert(0, xp)
    from psmgrating import twowave as tw
    PROBE = {"lumerical_R": (0.0, -14.0, 7.5, -28.0, 20.0), "zemax_kba_R": (-45.0, -22.5, 0.0, 10.0),
             "lumerical_T": (0.0, 5.0, 10.0), "unslanted_T": (33.8627180249, 36.0, 30.0),
             "unslanted_T_rot5": (38.8627180249, 41.0, 35.0)}
    for name, c in CASES.items():
        thr = float(np.arccos(c["K"] * c["lam_r_um"] * 1e3 / (4 * np.pi * c["n0"])))
        for th in PROBE[name]:
            k0 = 2 * np.pi / c["lam_nm"]
            psi = c["psi"] if np.cos(th * deg + c["psi"]) >= 0 else c["psi"] + np.pi
            thc = th * deg + psi
            for pol, xpol, chi0, chi1 in (("TE", "sigma", 3e-4, 1e-4), ("TM", "pi", 0.0, 0.0), ("TM", "pi", 3e-4, 1e-4),
                                          ("TE", "sigma", 1e-3, 0.0), ("TM", "pi", 0.0, 2e-4)):
                mine = float(kogelnik.efficiency(k0, c["n0"], c["n1"], c["K"], psi, thc, c["d_nm"], pol=pol, chi0=chi0, chi1=chi1))
                x = float(tw.efficiency("kogelnik", xpol, n0=c["n0"], n1=c["n1"], d=c["d_nm"] * 1e-9, lam_r=c["lam_r_um"] * 1e-6,
                                        theta_r=thr, lam_c=c["lam_nm"] * 1e-9, theta_c=thc, psi=psi, chi0=chi0, chi1=chi1))
                check(f"K3 {name} th={th:+.2f} {pol} chi0={chi0:g} chi1={chi1:g}", mine, x, 1e-10)
except ImportError as e:
    print("SKIP K3 (explorer package not importable -- expected on a machine without the psm-explorer checkout):", e)
    skips.append("K3")

print("K4  TM zero-coupling locus")
c = CASES["lumerical_R"]
k0 = 2 * np.pi / c["lam_nm"]
v = float(kogelnik.efficiency(k0, c["n0"], c["n1"], c["K"], c["psi"], 45 * deg, c["d_nm"], pol="TM"))
check("K4 TM eta at theta_c = 45 deg", v, 0.0, 1e-30)

print(f"\n{len(fails)} failure(s); {len(skips)} gate(s) skipped for want of an oracle on this machine: {skips}")
print("VERIFIED 2026-10-16 on the benchmark PC with the explorer present: 0 failures, 0 skips (evidence/kogelnik_gates.txt)")
sys.exit(1 if fails else 0)
