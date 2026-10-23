"""Check fpsm.py against every number in Section 6 of the F-PSM note.

Run:  python test_fpsm_vectors.py      (exit code 1 on any failure)
The note prints 10 significant figures; tolerance is 5e-10 relative on the
printed value unless the note itself marks the entry approximate.
"""
import sys

import numpy as np

import fpsm

deg = np.pi / 180.0
fails = []


def check(name, got, want, rtol=5e-10, atol=0.0):
    err = abs(got - want)
    ok = err <= atol + rtol * abs(want)
    print(f"{'PASS' if ok else 'FAIL'}  {name:34s} got {got!r:>28}  want {want!r}")
    if not ok:
        fails.append(name)


# ---- 6.1 transmission, dimensionless -------------------------------------
k0 = 1.0 / np.cos(30 * deg)
g_want = 2.3094010768e-2
for gd, d, eta in [(0.5, 21.6506350946, 0.2301656233),
                   (1.5, 64.9519052838, 0.9935456466),
                   (3.0, 129.9038105677, 0.0196892436)]:
    r = fpsm.fpsm_detail(k0, 1.0, 0.02, 2.0, 90 * deg, 30 * deg, d)
    check(f"T dimless gd={gd} case", int(r["case"]), 1)
    check(f"T dimless gd={gd} gamma0", float(r["gamma0"]), 5.7735026919e-1)
    check(f"T dimless gd={gd} g", float(r["g"]), g_want)
    check(f"T dimless gd={gd} ge", float(r["ge"]), 2.3100682371e-2)
    check(f"T dimless gd={gd} eta", float(r["eta"]), eta, rtol=0, atol=6e-11)   # printed to 10 d.p.

# ---- 6.1 transmission, nanometres ----------------------------------------
k0 = 2 * np.pi / 500.0
for tag, psi, thc, K, d, g0, g1, ge, eta in [
        ("T-a", 90, 123.8627180249, -2.1006130130e-2, 10000, 1.5652200570e-2, 1.5652200570e-2,
         2.2696851754e-4, 0.5860353852),
        ("T-b", 80, 113.8627180249, -1.5251047602e-2, 5000, 1.5652200570e-2, 1.8300517193e-2,
         2.0986734470e-4, 0.7518324125)]:
    r = fpsm.fpsm_detail(k0, 1.5, 0.03, K, psi * deg, thc * deg, d)
    check(f"{tag} case", int(r["case"]), 1)
    check(f"{tag} gamma0", float(r["gamma0"]), g0)
    check(f"{tag} gamma1", float(r["gamma1"]), g1)
    check(f"{tag} ge", float(r["ge"]), ge)
    check(f"{tag} eta", float(r["eta"]), eta, rtol=0, atol=2e-9)   # K printed to 11 s.f.
r = fpsm.fpsm_detail(k0, 1.5, 0.03, -1.5251047602e-2, 80 * deg, 113.8627180249 * deg, 5000)
check("T-b delta", float(r["delta"]), 5.1702028204e-7, rtol=1e-3)  # difference of ~1e-2 numbers

# ---- 6.2 reflection, dimensionless ---------------------------------------
k0 = 1.0 / np.cos(25 * deg)
for d, ReBm, ImBm, eta in [(39.3593963, -7.7818896e-6, 5.5801350e-7, 0.2134560937),
                           (118.0781889, -7.6226552e-6, 1.6626223e-6, 0.8190888015),
                           (236.1563778, 7.0932424e-6, -3.2488610e-6, 0.9899934123)]:
    r = fpsm.fpsm_detail(k0, 1.0, 0.02, 2.0, 15 * deg, 25 * deg, d)
    check(f"R dimless d={d:.1f} case", int(r["case"]), -1)
    check(f"R dimless d={d:.1f} a", float(r["gamma0"]), 1.0866151291)
    check(f"R dimless d={d:.1f} b", float(r["gamma1"]), 8.4523652348e-1)
    check(f"R dimless d={d:.1f} g", float(r["g"]), 1.2703446876e-2)
    check(f"R dimless d={d:.1f} ge", float(r["ge"]), 1.2701192970e-2)
    check(f"R dimless d={d:.1f} delta", float(r["delta"]), 8.2147551630e-5, rtol=1e-8)
    check(f"R dimless d={d:.1f} Re Hs", complex(r["Hs"]).real, 6.5988346166e-5, rtol=1e-9)
    check(f"R dimless d={d:.1f} Im Hs", complex(r["Hs"]).imag, 3.7920252675e-5, rtol=1e-9)
    check(f"R dimless d={d:.1f} Re B-", complex(r["Bm"]).real, ReBm, rtol=1e-6)   # d printed to 9 s.f.
    check(f"R dimless d={d:.1f} Im B-", complex(r["Bm"]).imag, ImBm, rtol=2e-5)
    check(f"R dimless d={d:.1f} eta", float(r["eta"]), eta, rtol=0, atol=5e-10)   # d printed to 10 s.f.: d(eta)/dd * 5e-8 ~ 3e-10

# ---- 6.2 reflection, nanometres ------------------------------------------
k0 = 2 * np.pi / 500.0
for tag, psi, thc, K, d, a, b, ge, delta, ReH, ImH, eta in [
        ("R-a", -15, 18.8627180249, 3.5674516079e-2, 10000, 1.5652200570e-2, 1.8806735851e-2,
         2.0705765133e-4, 1.4051464544e-6, 5.6014195772e-5, 5.2714659890e-5, 0.9382538133),
        ("R-b", -25, 8.8627180249, 3.7248996723e-2, 5000, 1.5652200570e-2, 1.8106855219e-2,
         2.1102129040e-4, 1.3315646491e-6, 3.9709993011e-5, 3.7688681491e-5, 0.6142648776)]:
    r = fpsm.fpsm_detail(k0, 1.5, 0.03, K, psi * deg, thc * deg, d)
    check(f"{tag} case", int(r["case"]), -1)
    check(f"{tag} a", float(r["gamma0"]), a)
    check(f"{tag} b", float(r["gamma1"]), b)
    check(f"{tag} ge", float(r["ge"]), ge)
    check(f"{tag} delta", float(r["delta"]), delta, rtol=1e-4)
    check(f"{tag} Re Hs", complex(r["Hs"]).real, ReH, rtol=1e-8)
    check(f"{tag} Im Hs", complex(r["Hs"]).imag, ImH, rtol=1e-8)
    check(f"{tag} eta", float(r["eta"]), eta, rtol=0, atol=2e-9)

# ---- 6.3 structural checks -----------------------------------------------
for tag, args in [("T", (2 * np.pi / 500, 1.5, 0.03, -2.1006130130e-2, 90 * deg, 123.8627180249 * deg, 1e4)),
                  ("R", (2 * np.pi / 500, 1.5, 0.03, 3.5674516079e-2, -15 * deg, 18.8627180249 * deg, 1e4))]:
    z = list(args); z[2] = 0.0
    check(f"{tag} n1=0 -> eta=0", float(fpsm.efficiency(*z)), 0.0, atol=0.0)
    z = list(args); z[6] = 0.0
    check(f"{tag} d=0 -> eta=0", float(fpsm.efficiency(*z)), 0.0, atol=1e-30)

# ---- 2.2 sign recovery -----------------------------------------------------
check("sign recovery T-a", fpsm.recover_K_sign(2 * np.pi / 500, 1.5, 2.1006130130e-2, 90 * deg,
                                               123.8627180249 * deg), -2.1006130130e-2)
check("sign recovery R-a", fpsm.recover_K_sign(2 * np.pi / 500, 1.5, 3.5674516079e-2, -15 * deg,
                                               18.8627180249 * deg), 3.5674516079e-2)

print(f"\n{len(fails)} failure(s)")
sys.exit(1 if fails else 0)
