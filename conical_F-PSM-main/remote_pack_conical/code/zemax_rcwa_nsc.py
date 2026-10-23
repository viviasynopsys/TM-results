"""Drive the Ansys Zemax RCWA volume-hologram DLL (``vhg_constant_hologram230313.dll``,
Knowledgebase article 42661816933139) on the same index-matched NSC bench that
``zemax_kogelnik_nsc.py`` uses, and tabulate first-order efficiency vs replay
angle and vs wavelength.

STATUS: **VERIFIED 2026-09-18** on Ansys Zemax OpticStudio 2024 R1.00, PremiumEdition seat, standalone
(headless) ZOS-API, DLL vhg_constant_hologram230313.dll in <ZemaxDataDir>\DLL\Diffractive.
Probe printouts: ../results/zemax_rcwa_probe.txt (gates Z1-Z4) and ../results/zemax_rcwa_probe2.txt
(order filter, polarisation, 160/320 layers, thickness).  What was settled, in short:
  * parameter names / 0-based indices / template values: exactly as the comment below (gate Z1);
  * "Holo Thic (um)" is micrometres (6 -> the 6 um response; 5.84/6.04/6.2 move eta as the
    referee's thickness derivative says) (gate Z2);
  * "Only these orders" = 16 (shipped) is an ORDER FILTER: order 0 and all transmitted orders read
    zero; with 0 the energy over orders -5..5 is 1.0000 and order 1 is unchanged (gate Z3);
  * Max Order 5 = 9 = 13 to 1e-10; Layer per period 10 -> 8e-3 low, 160 vs 320 <= 2.5e-5 (gate Z4);
  * "Interpolation" is an angular lookup-grid size (0 -> 201 points); 0 vs 16001 changes <= 3e-6;
  * the DLL is polarisation-aware (TM differs); Jx = 1 is TE on this y-z bench;
  * the DLL expands the permittivity as eps = n0^2 + 2 n0 dn cos(K.r), dropping the dn^2 terms,
    although its article documents n = n0 + dn cos(K.r): rcwa_te with eps_form="linear"
    reproduces the converged DLL to ~1e-5 where the exact form differs by 3.5e-3 (theta = -14 deg,
    lumerical_R).  Nothing here is tuned to hide that; the sweeps use dn = n1, Holo index = n0.

Method
    1. Read every Diffraction-tab parameter of the article's own, working file
       ``../zemax/ansys_kba/test_example1/test_example1.zmx`` (name -> value, both
       the Reflect and the Transmit column).  Unknown parameters (Test Mode, Only
       these orders, Stochastic mode, Layer per period, ...) therefore keep the
       values Ansys shipped -- we never guess them.
    2. Build the bench of ``zemax_kogelnik_nsc.build(volume_hologram=False)``.
    3. Point the Hologram Lens' Diffraction tab at the DLL, write the template
       values back, then override BY NAME only what defines the grating.
    4. Start Order = Stop Order = the first order, Split NSC Rays ON: the entry-side
       detector then reads R(+1), the far-side detector T(+1), in watts per 1 W.

Usage
    python zemax_rcwa_nsc.py --probe lumerical_R    # settle the assumptions, write nothing
    python zemax_rcwa_nsc.py lumerical_R            # -> ../results/zemax_rcwa_lumerical_R.csv
    python zemax_rcwa_nsc.py --pol TM --probe lumerical_R   # TM probe (Jy = 1): ladders at 10..320 layers
    python zemax_rcwa_nsc.py --pol TM lumerical_R   # -> ../results/zemax_rcwa_lumerical_R_TM.csv

--pol TM (added 2026-10-16; the switch itself VERIFIED 2026-09-18 by probe 2b, which gave
0.7272546 / 0.0245981 / 0.0520480 at theta = 0, -14, 7.5 deg on lumerical_R at 10 layers per
period against 0.7965196 / 0.0245959 / 0.0978717 for TE; the sweeps are UNTESTED).  Output file
names carry the suffix _TM and the saved bench is zemax/rcwa_dll_<case>_TM.zmx.  The probe's
layer ladder now runs 10..320 so the TM convergence is on record like the TE one.
"""
import csv
import math
import os
import re
import sys

from cases import CASES
from zemax_kogelnik_nsc import ROOT, build, trace
from zos_connect import Zos

DLL = "vhg_constant_hologram230313.dll"
TEMPLATE = os.path.join(ROOT, "zemax", "ansys_kba", "test_example1", "test_example1.zmx")

# DLL parameter names, in the DLL's own order (strings read out of the binary):
#   Max Order | X1 Y1 Z1 X2 Y2 Z2 | Const. Wave (um) | Holo Thic (um) | n1 | n2 | Holo index |
#   Index modulation (dn) | Use Coat File | Shrinkage (unused) | Index Shift (unused) |
#   Rotate Grating | Interpolation | Test Mode | Only these orders | Stochastic mode |
#   Layer per period | Unused | Holo Type
# Values in test_example1.zmx, same order (the file's NSOF 101..124 records):
#   5 | 0 1 -1 0 0 -1 | 0.55 | 5 | 1.5 | 1.5 | 1.5 | 0.05 | 0 | 0 | 0 | 0 | 0 | 1 | 16 | 0 | 10 | 0 | 2
# ASSUMED: that pairing of names to values (the binary lists 24 names, the file 24 non-trivial
# values, but the API is the authority -- read the names live, never trust this comment).


def diffraction_params(dd):
    """{name: (index, reflect_value, transmit_value)} read live.  Indices are 0-based
    (VERIFIED 2026-09-17 with Diff2DSample.dll: names at 0..N-1, '' at N); the scan to N
    and the empty-name filter stay as a guard.  Note the API's own spelling,
    ``GetTransmitParamaterName``."""
    out = {}
    for i in range(0, int(dd.NumberOfParameters) + 1):
        try:
            name = str(dd.GetReflectParameterName(i)).strip()
        except Exception:
            continue
        if name and name not in out:
            out[name] = (i, float(dd.GetReflectParameterValue(i)), float(dd.GetTransmitParameterValue(i)))
    return out


def read_template(z):
    if not os.path.isfile(TEMPLATE):
        raise FileNotFoundError(f"{TEMPLATE} missing -- restore test_example1.zar first (see execution.md)")
    if not z.sys.LoadFile(TEMPLATE, False):
        raise RuntimeError("OpticStudio refused to load the article's test_example1.zmx")
    holo = z.sys.NCE.GetObjectAt(2)
    dd = holo.DiffractionData
    info = dict(split=str(dd.Split), dll=str(dd.DLL), start=int(dd.StartOrder), stop=int(dd.StopOrder),
                params=diffraction_params(dd))
    if not info["params"]:
        raise RuntimeError("template exposes no diffraction parameters -- DLL not loaded?")
    return info


def find(params, *needles):
    """Parameter name containing every needle (case-insensitive); loud if not unique."""
    hits = [n for n in params if all(s.lower() in n.lower() for s in needles)]
    if len(hits) != 1:
        raise KeyError(f"expected exactly one parameter matching {needles}, got {hits}; have {list(params)}")
    return hits[0]


# Sweep settings.  Article 42661761082899 (the sibling dynamic-link DLL, same "Interpolation"
# semantics -- this DLL's own strings say "interp_samp" and "Rays' incident angle outside of
# sampling grids"): Interpolation is the size of an angular LOOKUP GRID, not an on/off switch.
# 0 or 1 -> 201 samples; 2..50 -> 51; > 51 -> that many samples; negative -> the same grid
# with power-conserving interpolation.  The old "0 = exact per-ray solve" reading was wrong.
# The defaults below are what --probe converged on (see results/zemax_rcwa_probe.txt); the
# CLI flags --interp, --layers, --max-order override them and are recorded in the CSV.
INTERP_DEFAULT = 0.0        # the 201-point grid.  16001 gave the same values (<= 3e-6 at the probe angles, and
                            # <= the value in results/zemax_rcwa_verdict.md over the whole lumerical_R angle sweep,
                            # results/zemax_rcwa_lumerical_R.angle_only_run1.csv) but took 36 s per first call and
                            # KILLED the OpticStudio server at every wavelength change (results/zemax_rcwa_probe3.txt).
LAYERS_DEFAULT = 160.0      # z-slices per period: 160 vs 320 moves eta by <= 2.5e-5 (gate Z4 asks < 1e-4); 10 -> 8e-3 off
MAXORDER_DEFAULT = 5.0      # harmonics: 5 vs 9 vs 13 identical to 1e-10 (gate Z4); the template's own value
# VERIFIED 2026-09-18 (results/zemax_rcwa_probe2.txt): the template's "Only these orders" = 16 is an
# order FILTER -- with it, order 0 and every transmitted order return zero power (energy sum over
# orders -5..5 = 0.797).  With 0 the sum is 1.0000 and the order-1 value is unchanged, i.e. no
# renormalisation.  "Test Mode" (template 1) had no effect on any traced value; kept at 1.
DLL_OVERRIDES = {"Only these orders": 0.0}


def configure(z, holo, case, template, order, max_order=None, layers=None, interp=None):
    Z = z.ZOSAPI
    dd = holo.DiffractionData
    dd.Split = Z.Editors.NCE.DiffractionSplitType.SplitByDLL
    dd.DLL = DLL
    live = diffraction_params(dd)
    if not live:
        raise RuntimeError(f"{DLL} exposes no parameters -- is it in <Zemax data>\\DLL\\Diffractive ?")
    for name, (_, r, t) in template["params"].items():           # Ansys' shipped values first
        if name in live:
            dd.SetReflectParameterValue(live[name][0], r)
            dd.SetTransmitParameterValue(live[name][0], t)
    far = max(abs(v) for v in case["P1"] + case["P2"])
    grating = {                                                   # then the grating, by name
        find(live, "X1"): case["P1"][0] / far, find(live, "Y1"): case["P1"][1] / far,
        find(live, "Z1"): case["P1"][2] / far, find(live, "X2"): case["P2"][0] / far,
        find(live, "Y2"): case["P2"][1] / far, find(live, "Z2"): case["P2"][2] / far,
        find(live, "Wave"): case["lam_r_um"],
        find(live, "Thic"): case["d_um"],                         # MICROMETRES: settled at gate Z2 (see probe printout)
        find(live, "n1"): case["n0"], find(live, "n2"): case["n0"],
        find(live, "Holo index"): case["n0"],
        find(live, "modulation"): case["n1"],
        find(live, "Holo Type"): 2.0,
        find(live, "Interpolation"): float(INTERP_DEFAULT if interp is None else interp),
    }
    if max_order is not None:
        grating[find(live, "Max Order")] = float(max_order)
    if layers is not None:
        grating[find(live, "Layer per period")] = float(layers)
    for k, v in DLL_OVERRIDES.items():
        grating[find(live, k)] = float(v)
    for name, v in grating.items():
        dd.SetReflectParameterValue(live[name][0], float(v))
        dd.SetTransmitParameterValue(live[name][0], float(v))
    dd.StartOrder = int(order)
    dd.StopOrder = int(order)
    return diffraction_params(dd)


def flag(argv, key, default, cast=float):
    """--key VALUE  (absent -> default)."""
    if key in argv:
        return cast(argv[argv.index(key) + 1])
    return default


def dll_log():
    """The DLL writes vhg_constant_hologram.txt beside itself -- return its tail, or why not."""
    for d in (os.path.join(os.path.expanduser("~"), "Documents", "Zemax", "DLL", "Diffractive"),):
        p = os.path.join(d, "vhg_constant_hologram.txt")
        if os.path.isfile(p):
            with open(p, errors="replace") as f:
                lines = f.readlines()
            return f"{p}: {len(lines)} lines\n" + "".join(lines[-40:])
    return "(no vhg_constant_hologram.txt found beside the DLL)"


def main(argv):
    probe = "--probe" in argv
    pol = "TE"
    if "--pol" in argv:                                   # parse (and remove) before the case name is looked for
        i = argv.index("--pol")
        pol = argv[i + 1].upper()
        del argv[i:i + 2]
    if pol not in ("TE", "TM"):
        raise SystemExit(f"--pol must be TE or TM, not {pol!r}")
    suf = "" if pol == "TE" else "_TM"
    name = [a for a in argv if not a.startswith("--") and not re.fullmatch(r"-?\d+(\.\d+)?", a)][0]
    case = CASES[name]
    interp = flag(argv, "--interp", INTERP_DEFAULT)
    layers = flag(argv, "--layers", LAYERS_DEFAULT)
    max_order = flag(argv, "--max-order", MAXORDER_DEFAULT)
    z = Zos()
    print(f"OpticStudio {z.app.ZOSMajorVersion}.{z.app.ZOSMinorVersion}.{z.app.ZOSSPVersion}  licence {z.app.LicenseStatus}  data {z.app.ZemaxDataDir}  pol {pol}", flush=True)
    dll_path = os.path.join(str(z.app.ZemaxDataDir), "DLL", "Diffractive", DLL)
    print(f"DLL present at {dll_path}: {os.path.isfile(dll_path)}", flush=True)
    rows, failed = [], False
    kw = dict(interp=interp, layers=layers, max_order=max_order)
    try:
        template = read_template(z)
        print("template:", {k: v for k, v in template.items() if k != "params"}, flush=True)
        for n, (i, r, t) in template["params"].items():
            print(f"   [{i:2d}] {n:28s} reflect={r:g} transmit={t:g}", flush=True)

        nce, holo, src = build(z, case, volume_hologram=False, pol=pol)
        w1 = z.sys.SystemData.Wavelengths.GetWavelength(1)
        src.TiltAboutX = float(case["zemax_tilt_bragg_deg"])

        # Which integer is "the" first order?  Decided by measurement at Bragg (both signs traced).
        at_bragg = {}
        for order in (1, -1):
            live = configure(z, holo, case, template, order, **kw)
            at_bragg[order] = trace(z, nce, split=True)
        print("live parameters after configure():", flush=True)
        for n, (i, r, t) in live.items():
            print(f"   [{i:2d}] {n:28s} reflect={r:g} transmit={t:g}", flush=True)
        print("at Bragg (entry-side W, far-side W) by order:", at_bragg, flush=True)
        order = max(at_bragg, key=lambda k: max(at_bragg[k]))
        if max(at_bragg[order]) < 0.3:
            raise RuntimeError(f"no order carries the Bragg peak: {at_bragg} -- geometry or parameter mapping is wrong")

        if probe:
            total = [0.0, 0.0]
            for m in range(-3, 4):                                # energy book-keeping across orders (gate Z3)
                configure(z, holo, case, template, m, **kw)
                pr, pt = trace(z, nce, split=True)
                total[0] += pr
                total[1] += pt
                print(f"   order {m:+d}: R={pr:.6f} T={pt:.6f}", flush=True)
            print(f"   sum over orders -3..3: R={total[0]:.6f} T={total[1]:.6f} total={sum(total):.6f}  (expect 1.000 +- 1e-3)", flush=True)
            # Probe angles: Bragg plus the case's off-Bragg reference angles of execution.md section 6.
            off = {"lumerical_R": (-14.0, 7.5), "lumerical_T": (5.0,), "zemax_kba_R": (-22.5,), "unslanted_T": (36.0,)}[name]
            taus = [case["zemax_tilt_bragg_deg"]] + [-t for t in off]
            def at(taus_, **k):
                out = []
                for tau in taus_:
                    src.TiltAboutX = float(tau)
                    out.append(trace(z, nce, split=True))
                src.TiltAboutX = float(case["zemax_tilt_bragg_deg"])
                return out
            print(f"   probe angles theta = {[-t for t in taus]} deg  (tilt = {taus})", flush=True)
            for mo in (3, 5, 9, 13):                              # harmonic convergence (gate Z4)
                configure(z, holo, case, template, order, interp=interp, layers=layers, max_order=mo)
                print(f"   Max Order {mo:2d}: (R, T) = {at(taus)}", flush=True)
            for lp in (10, 20, 40, 80, 160, 320):                 # z-slicing convergence (gate Z4, second half; 160/320 added 2026-10-16)
                configure(z, holo, case, template, order, interp=interp, layers=lp, max_order=max_order)
                print(f"   Layer per period {lp:3d}: (R, T) = {at(taus)}", flush=True)
            for ip in (0, 1001, 4001, 16001, -4001):              # angular lookup-grid convergence
                configure(z, holo, case, template, order, interp=ip, layers=layers, max_order=max_order)
                print(f"   Interpolation {ip:6d}: (R, T) = {at(taus)}", flush=True)
            print("DLL log tail:\n" + dll_log(), flush=True)
            return

        configure(z, holo, case, template, order, **kw)
        print(f"sweep settings: order={order} interp={interp} layers={layers} max_order={max_order} pol={pol}", flush=True)
        os.makedirs(os.path.join(ROOT, "zemax"), exist_ok=True)
        z.sys.SaveAs(os.path.join(ROOT, "zemax", f"rcwa_dll_{name}{suf}.zmx"))
        # Order labelling.  The DLL's integer orders are fixed; OpticStudio's native model (and
        # compare.py, for the Python models) calls "order 1" the order with K.k_playback > 0, i.e.
        # cos(theta + psi) >= 0, which swaps to the conjugate past the fringe-parallel angle
        # (lumerical_T, theta > 14 deg).  Columns P_entry/P_far follow that rule so the table is
        # comparable; the raw values of BOTH integer orders are kept (P_*_plus / P_*_minus) wherever
        # the rule flips inside the sweep, so nothing is hidden.
        psi = float(case["psi"])
        flips = [math.cos(-tau * math.pi / 180.0 + psi) < 0 for tau in case["zemax_tilts_deg"]]
        both = any(flips)
        print(f"order rule: {sum(flips)} of {len(flips)} angle points use the conjugate order {-order}; both orders traced: {both}", flush=True)
        for tau, flip in zip(case["zemax_tilts_deg"], flips):
            src.TiltAboutX = float(tau)
            w1.Wavelength = case["lam_um"]
            if both:
                val = {}
                for o in (order, -order):
                    dd = holo.DiffractionData
                    dd.StartOrder = int(o)
                    dd.StopOrder = int(o)
                    val[o] = trace(z, nce, split=True)
                used = -order if flip else order
                pr, pt = val[used]
                plus, minus = val[abs(order)], val[-abs(order)]
                rows.append(("angle", tau, case["lam_um"], pr, pt, used, interp, layers, max_order,
                             plus[0], plus[1], minus[0], minus[1]))
            else:
                pr, pt = trace(z, nce, split=True)
                rows.append(("angle", tau, case["lam_um"], pr, pt, order, interp, layers, max_order, "", "", "", ""))
        if both:
            dd = holo.DiffractionData
            dd.StartOrder = int(order)
            dd.StopOrder = int(order)
        src.TiltAboutX = float(case["zemax_tilt_bragg_deg"])
        for lam in case["zemax_lams_um"]:
            w1.Wavelength = float(lam)
            pr, pt = trace(z, nce, split=True)
            rows.append(("wavelength", case["zemax_tilt_bragg_deg"], lam, pr, pt, order, interp, layers, max_order, "", "", "", ""))
        w1.Wavelength = case["lam_um"]
        z.sys.Save()
        print("DLL log tail:\n" + dll_log(), flush=True)
    except Exception:
        import traceback
        traceback.print_exc()
        failed = True
    finally:
        if not probe:
            bad = failed or not rows
            os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
            out = os.path.join(ROOT, "results", f"zemax_rcwa_{name}{suf}{'.FAILED' if bad else ''}.csv")
            with open(out, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["sweep", "source_tilt_x_deg", "lambda_um", "P_entry_side_W", "P_far_side_W", "order",
                            "interpolation", "layer_per_period", "max_order",
                            "P_entry_plus_W", "P_far_plus_W", "P_entry_minus_W", "P_far_minus_W"])
                w.writerows(rows)
            print(f"{len(rows)} rows -> {out}", flush=True)
        try:
            z.close()
        except Exception as e:
            print("close:", e)
        sys.stdout.flush()
        os._exit(1 if failed else 0)       # pythonnet + ZOS-API hangs on normal interpreter exit


if __name__ == "__main__":
    main(sys.argv[1:])
