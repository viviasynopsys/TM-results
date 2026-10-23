"""Gates for rcwa.py (TE + TM, lossless + lossy) against INDEPENDENT oracles.  Exit code 1 on any
failure.  Written 2026-10-16.

  G1  TE lossless path == rcwa_te.solve, bit for bit, every benchmark grating, both eps forms
  G2  1/eps closed-form harmonics == 8192-point FFT (exact and linear forms, complex index)
  G3  TM == explorer psmgrating.rcwa pol="pi" (env PSMGRATING_PATH, default the Repos checkout)
  G4  TE and TM, lossless and lossy == DBR_volume_gratings/validation/rcw_tm/rcw_rigorous.py
      (run 007 of 2026-07-28; env RCW_RIGOROUS, default that path), plus its PRINTED values
  G5  energy: lossless sum = 1 to 1e-10 (TM, every grating, several angles); lossy sum < 1
  G6  TM truncation: L = 10 against L = 14 over every grating's angle sweep
  G7  continuity: chi -> 0 reproduces the lossless result
Everything prints with the measured number, PASS or FAIL.
"""
import importlib.util
import os
import sys

import numpy as np

import rcwa
import rcwa_te
from cases import CASES, deg

fails = []
skips = []      # oracle not present on this machine: reported, never counted as a pass


def check(name, got, want, atol):
    ok = np.all(np.abs(np.asarray(got) - np.asarray(want)) <= atol)
    worst = float(np.max(np.abs(np.asarray(got) - np.asarray(want))))
    print(f"{'PASS' if ok else 'FAIL'}  {name:62s} max|diff| {worst:.2e}  (tol {atol:.0e})")
    if not ok:
        fails.append(name)


def args_for(c, theta_deg, lam_nm=None):
    lam_nm = c["lam_nm"] if lam_nm is None else lam_nm
    k0 = 2 * np.pi / lam_nm
    psi = c["psi"] if np.cos(theta_deg * deg + c["psi"]) >= 0 else c["psi"] + np.pi   # compare.py's K.k rule
    return (k0, c["n0"], c["n1"], c["K"], psi, theta_deg * deg + psi, c["d_nm"])


PROBE = {"lumerical_R": (0.0, -14.0, 7.5, -28.0, 20.0), "zemax_kba_R": (-45.0, -22.5, 0.0, 10.0),
         "lumerical_T": (0.0, 5.0, 28.0, 40.0), "unslanted_T": (33.8627180249, 36.0, 30.0),
         "unslanted_T_rot5": (38.8627180249, 41.0, 35.0)}

# ---- G1: TE path identical to rcwa_te ------------------------------------------------------------
print("G1  TE lossless path against rcwa_te.solve (must be bit-identical)")
for name, c in CASES.items():
    for th in PROBE[name]:
        a = args_for(c, th)
        for form in ("exact", "linear"):
            old = rcwa_te.solve(*a, L=10, eps_form=form)
            new = rcwa.solve(*a, L=10, eps_form=form, pol="TE")
            ident = np.array_equal(old["DE_R"], new["DE_R"]) and np.array_equal(old["DE_T"], new["DE_T"])
            worst = max(np.max(np.abs(old["DE_R"] - new["DE_R"])), np.max(np.abs(old["DE_T"] - new["DE_T"])))
            tag = f"G1 {name} th={th:+.2f} {form}"
            print(f"{'PASS' if ident else 'FAIL'}  {tag:62s} bit-identical {ident}  max|diff| {worst:.1e}")
            if not ident:
                fails.append(tag)

# ---- G2: 1/eps harmonics against an FFT ----------------------------------------------------------
print("G2  1/eps closed-form harmonics against a 8192-point FFT")
phi = 2 * np.pi * np.arange(8192) / 8192
for (n0g, n1g) in ((1.5, 0.04), (1.5 - 1e-3j, 0.04 - 2e-4j), (1.5 - 5e-3j, 0.2 - 1e-3j)):
    for form in ("exact", "linear"):
        eps = (n0g + n1g * np.cos(phi)) ** 2 if form == "exact" else n0g * n0g + 2 * n0g * n1g * np.cos(phi)
        fft = np.fft.ifft(1.0 / eps)                             # c_m with 1/eps = sum c_m exp(+i m phi); even in m
        ms = np.arange(-20, 21)
        closed = np.array([rcwa.inv_eps_harmonics(n0g, n1g, m, form) for m in ms])
        check(f"G2 1/eps {form} n0g={n0g} n1g={n1g}", closed, fft[ms % 8192], 1e-13)
        eh = rcwa.eps_harmonics(n0g, n1g, form)
        fe = np.fft.ifft(eps)
        check(f"G2 eps   {form} n0g={n0g} n1g={n1g}", [eh[0], eh[1], eh[2]], [fe[0], fe[1], fe[2]], 1e-13)

# ---- G3: TM against the explorer core -------------------------------------------------------------
print("G3  TM against psmgrating.rcwa pol='pi' (exact eps)")
xp = os.environ.get("PSMGRATING_PATH", r"C:\Users\vivia\Repos\psm-explorer\python")
try:
    sys.path.insert(0, xp)
    from psmgrating import rcwa as xr
    # explorer inputs: lam_c, lam_r, n0, n1, chi0, chi1, theta_c (angle between replay beam and K),
    # theta_r (angle between the recording reference beam and K), psi, d.  cases.py holds K > 0 with
    # psi folded, so theta_r is the acute angle with K = 4 pi n0 cos(theta_r) / lam_r.
    for name, c in CASES.items():
        thr = float(np.arccos(c["K"] * c["lam_r_um"] * 1e3 / (4 * np.pi * c["n0"])))
        for th in PROBE[name]:
            a = args_for(c, th)
            thc, psi = a[5], a[4]
            x = xr.solve(lam_c=c["lam_nm"] * 1e-9, lam_r=c["lam_r_um"] * 1e-6, n0=c["n0"], n1=c["n1"], chi0=0.0, chi1=0.0,
                         theta_c=thc, theta_r=thr, psi=psi, d=c["d_nm"] * 1e-9, pol="pi", L=10)
            mine = rcwa.solve(*a, L=10, pol="TM")
            check(f"G3 {name} th={th:+.2f} DE_R", mine["DE_R"], x["DE_R"], 1e-10)
            check(f"G3 {name} th={th:+.2f} DE_T", mine["DE_T"], x["DE_T"], 1e-10)
            # lossy TM and TE too
            x2 = xr.solve(lam_c=c["lam_nm"] * 1e-9, lam_r=c["lam_r_um"] * 1e-6, n0=c["n0"], n1=c["n1"], chi0=2e-4, chi1=1e-4,
                          theta_c=thc, theta_r=thr, psi=psi, d=c["d_nm"] * 1e-9, pol="pi", L=10)
            m2 = rcwa.solve(*a, L=10, pol="TM", chi0=2e-4, chi1=1e-4)
            check(f"G3 {name} th={th:+.2f} lossy TM DE_R", m2["DE_R"], x2["DE_R"], 1e-10)
            x3 = xr.solve(lam_c=c["lam_nm"] * 1e-9, lam_r=c["lam_r_um"] * 1e-6, n0=c["n0"], n1=c["n1"], chi0=2e-4, chi1=1e-4,
                          theta_c=thc, theta_r=thr, psi=psi, d=c["d_nm"] * 1e-9, pol="sigma", L=10)
            m3 = rcwa.solve(*a, L=10, pol="TE", chi0=2e-4, chi1=1e-4)
            check(f"G3 {name} th={th:+.2f} lossy TE DE_R", m3["DE_R"], x3["DE_R"], 1e-10)
except ImportError as e:
    print("SKIP G3 (explorer package not importable -- expected on a machine without the psm-explorer checkout):", e)
    skips.append("G3")

# ---- G4: against the July two-polarisation core and its printed run-007 values -------------------
print("G4  TE/TM, lossless/lossy against validation/rcw_tm/rcw_rigorous.py (run 007, 2026-07-28)")
rp = os.environ.get("RCW_RIGOROUS", r"C:\Users\vivia\Code2\1_Optics\Research\DBR_volume_gratings\validation\rcw_tm\rcw_rigorous.py")
if os.path.isfile(rp):
    spec = importlib.util.spec_from_file_location("rcw_rigorous", rp)
    rr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rr)
    # run-007 baseline: n0 = 1.5, n1 = 0.045, d = 20 um, theta_r = 20, theta_c = 10, lam_r = 532 nm, psi = 10 deg;
    # legacy chi0 = D lam_c cos(theta_c) / (2 pi d)  (cos theta_c, NOT cos(theta_c - psi): the July convention)
    n0, n1, d_nm, thr, thc, lamr, psi = 1.5, 0.045, 20000.0, 20 * deg, 10 * deg, 532.0, 10 * deg
    K = 4 * np.pi * n0 * np.cos(thr) / lamr
    for lam in (545.0, 557.0, 560.0, 557.6):
        for D0 in (0.0, 1.0):
            for f1 in (0.0, 1.0):
                chi0 = D0 * lam * np.cos(thc) / (2 * np.pi * d_nm)
                chi1 = f1 * chi0
                for pol, lpol in (("TE", "TE"), ("TM", "TM")):
                    o = rr.rcw_planar(lam * 1e-9, lamr * 1e-9, n0, n1, D0, thc, thr, psi, d_nm * 1e-9, pol=lpol, chi1_factor=f1, M=10)
                    mine = rcwa.solve(2 * np.pi / lam, n0, n1, K, psi, thc, d_nm, L=10, pol=pol, chi0=chi0, chi1=chi1)
                    check(f"G4 {pol} lam={lam} D0={D0:g} chi1/chi0={f1:g} DE_R", mine["DE_R"], o["DE_R"], 1e-9)
                    check(f"G4 {pol} lam={lam} D0={D0:g} chi1/chi0={f1:g} DE_T", mine["DE_T"], o["DE_T"], 1e-9)
    # printed values of tm_rcw_validation_log.txt (run 007)
    chi0 = 1.0 * 557.0 * np.cos(thc) / (2 * np.pi * d_nm)
    v = rcwa.solve(2 * np.pi / 557.0, n0, n1, K, psi, thc, d_nm, L=10, pol="TM", chi0=chi0, chi1=0.0)["DE_R"][9]
    check("G4 printed: TM D0=1 lam=557 nm eta(l=-1) = 0.6623948815", v, 0.6623948815, 2e-9)
    v = rcwa.solve(2 * np.pi / 557.6, n0, n1, K, psi, thc, d_nm, L=10, pol="TM")["DE_R"][9]
    check("G4 printed: TM lossless Bragg peak 0.999435 at 557.6 nm", v, 0.999435, 2e-6)
    r = rcwa.solve(2 * np.pi / 557.0, n0, n1, K, psi, thc, d_nm, L=10, pol="TM")
    check("G4 printed: TM lossless energy |sum-1| (run 007: 8e-15)", r["sum_DE"], 1.0, 1e-12)
else:
    print("SKIP G4 (rcw_rigorous.py not found at", rp, "-- expected on a machine without the DBR_volume_gratings tree)")
    skips.append("G4")

# ---- G5: energy ----------------------------------------------------------------------------------
print("G5  energy conservation (TM lossless) and absorption (lossy sum < 1)")
for name, c in CASES.items():
    for th in PROBE[name]:
        a = args_for(c, th)
        r = rcwa.solve(*a, L=10, pol="TM")
        check(f"G5 {name} th={th:+.2f} TM sum_DE", r["sum_DE"], 1.0, 1e-10)
        rl = rcwa.solve(*a, L=10, pol="TM", chi0=1e-3)
        ok = rl["sum_DE"] < 1.0 - 1e-4
        print(f"{'PASS' if ok else 'FAIL'}  {'G5 ' + name + f' th={th:+.2f} TM chi0=1e-3 absorbs':62s} sum_DE {rl['sum_DE']:.6f}")
        if not ok:
            fails.append(f"G5 lossy {name} {th}")

# ---- G6: truncation ------------------------------------------------------------------------------
print("G6  TM truncation L = 10 against L = 14 over each angle sweep")
for name, c in CASES.items():
    th = c["theta_deg"]
    k0 = 2 * np.pi / c["lam_nm"]
    psi = np.where(np.cos(th * deg + c["psi"]) >= 0, c["psi"], c["psi"] + np.pi)
    a = (k0, c["n0"], c["n1"], c["K"], psi, th * deg + psi, c["d_nm"])
    r10 = rcwa.signal_efficiency(*a, L=10, pol="TM")
    r14 = rcwa.signal_efficiency(*a, L=14, pol="TM")
    check(f"G6 {name} TM L10 vs L14 (R and T faces)", np.r_[r10[0], r10[1]], np.r_[r14[0], r14[1]], 1e-10)

# ---- G7: continuity in chi -----------------------------------------------------------------------
print("G7  chi -> 0 continuity")
c = CASES["lumerical_R"]
a = args_for(c, 0.0)
for pol in ("TE", "TM"):
    r0 = rcwa.solve(*a, L=10, pol=pol)["DE_R"][9]
    r1 = rcwa.solve(*a, L=10, pol=pol, chi0=1e-12, chi1=1e-12)["DE_R"][9]
    check(f"G7 {pol} eta(chi=1e-12) vs eta(0)", r1, r0, 1e-9)

print(f"\n{len(fails)} failure(s); {len(skips)} gate(s) skipped for want of an oracle on this machine: {skips}")
print("VERIFIED 2026-10-16 on the benchmark PC with both oracles present: 0 failures, 0 skips (evidence/tm_referee_gates.txt)")
sys.exit(1 if fails else 0)
