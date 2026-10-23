"""Drive OpticStudio's NATIVE volume-hologram model (Kogelnik) in Non-Sequential
mode through ZOS-API and tabulate first-order efficiency vs replay angle and
vs wavelength, for the index-matched TE benches in ``cases.py``.

Bench (all inside one Rectangular Volume of model index n0, so there is no
Fresnel interface anywhere -- F-PSM's and the Lumerical article's assumption):

    1 Rectangular Volume   model glass n0
    2 Hologram Lens        same glass; Holo type 2; Volume Hologram? = 1;
                           Consider Fresnel? = 0; Diff Order = 1
    3 Source Ray           1 W, TE (Jx = 1: plane of incidence is y-z), tilt about X
    4 Detector Rectangle   ABSORB, entry side  (catches a reflected first order)
    5 Detector Rectangle   ABSORB, far side    (catches a transmitted first order)

Split NSC Rays is OFF and Diff Order = 1, so the only ray leaving face 1 is the
first order, carrying eta x 1 W; polarisation ON (the model needs it).

Usage:  python zemax_kogelnik_nsc.py <case-name> [--pol TE|TM]   -> ../results/zemax_<case>[_TM].csv
The .zmx is written by OpticStudio itself (never hand-edited) to ../zemax/.

--pol TM (added 2026-10-16, UNTESTED at OpticStudio): Source Ray Jx = 0, Jy = 1, i.e. the
electric field in the y-z plane of incidence.  Whether OpticStudio's native Kogelnik model
applies a polarisation factor at all is the question this run answers: compare the table with
kogelnik.efficiency(..., pol="TM") (kappa cos 2 theta_c) AND with pol="TE"; it will equal one
of them.  Both outcomes are findings; neither is a failure.
"""
import csv
import os
import sys

import numpy as np

from cases import CASES
from zos_connect import Zos

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def cell(o, header):
    for c in range(1, 80):
        ce = o.GetCellAt(c)
        if ce is None:
            break
        if str(ce.Header).strip() == header:
            return ce
    raise KeyError(f"no NCE column named {header!r} on {o.TypeName}")


def setd(o, header, v):
    """Set an NCE cell by column header.  The editor decides whether a column
    is Integer or Double (``Diff Order`` is Double, ``Holo type`` Integer), so
    ask the cell rather than guess."""
    ce = cell(o, header)
    if str(ce.DataType) == "Integer":
        ce.IntegerValue = int(round(v))
    else:
        ce.DoubleValue = float(v)


seti = setd


def set_model_glass(z, o, n0):
    mc = o.MaterialCell
    s = mc.CreateSolveType(z.ZOSAPI.Editors.SolveType.MaterialModel)
    s._S_MaterialModel.IndexNd = float(n0)
    s._S_MaterialModel.AbbeVd = 0.0
    s._S_MaterialModel.dPgF = 0.0
    mc.SetSolveData(s)


def build(z, case, volume_hologram=True, pol="TE"):
    """volume_hologram=False leaves the native Kogelnik model OFF (the RCWA DLL bench,
    ``zemax_rcwa_nsc.py``, supplies the efficiency through the Diffraction tab instead).
    pol: "TE" -> Jx = 1 (field along x, perpendicular to the y-z plane of incidence);
         "TM" -> Jy = 1 (field in the plane of incidence).  VERIFIED 2026-09-18 that the RCWA DLL
         responds to this switch (probe 2b); the native model's response is UNTESTED."""
    if pol not in ("TE", "TM"):
        raise ValueError(f"pol must be TE or TM, not {pol!r}")
    Z = z.ZOSAPI
    T = Z.Editors.NCE.ObjectType
    z.sys.New(False)
    z.sys.MakeNonSequential()
    z.sys.SystemData.Wavelengths.GetWavelength(1).Wavelength = case["lam_um"]
    # Default 1e-3 terminates any first-order ray with eta < 1e-3, which reads back
    # as a false zero in the sidelobes (VERIFIED 2026-09-17 on lumerical_R).
    z.sys.SystemData.NonSequentialData.MinimumRelativeRayIntensity = 1.0e-12
    nce = z.sys.NCE
    while nce.NumberOfObjects < 5:
        nce.AddObject()

    def make(i, typ):
        o = nce.GetObjectAt(i)
        o.ChangeType(o.GetObjectTypeSettings(typ))
        return o

    vol = make(1, T.RectangularVolume)
    vol.ZPosition = -2.0
    for h, v in [("X1 Half Width", 20), ("Y1 Half Width", 20), ("Z Length", 4.5),
                 ("X2 Half Width", 20), ("Y2 Half Width", 20)]:
        setd(vol, h, v)
    set_model_glass(z, vol, case["n0"])

    holo = make(2, T.HologramLens)
    set_model_glass(z, holo, case["n0"])
    setd(holo, "Radial Height", 5.0)
    setd(holo, "Thickness", 0.1)
    seti(holo, "Holo type", 2)
    seti(holo, "Diff Order", 1)
    setd(holo, "Cons. Wave", case["lam_r_um"])
    for h, v in zip(("X1", "Y1", "Z1", "X2", "Y2", "Z2"), case["P1"] + case["P2"]):
        setd(holo, h, v)
    if volume_hologram:
        seti(holo, "Volume Hologram?", 1)
        # LENS UNITS (mm), not micrometres: VERIFIED 2026-09-17 -- entering 6 gave the
        # Kogelnik response of a 6 mm plate to 6 digits at every sampled angle.
        setd(holo, "Hologram Thickness", case["d_um"] * 1.0e-3)
        setd(holo, "n1", case["n0"])
        setd(holo, "n2", case["n0"])
        setd(holo, "n", case["n0"])
        setd(holo, "dn", case["n1"])
        setd(holo, "Shrinkage", 0.0)
        setd(holo, "Index Shift", 0.0)
        seti(holo, "Consider Fresnel?", 0)
    else:
        seti(holo, "Volume Hologram?", 0)

    src = make(3, T.SourceRay)
    src.ZPosition = -0.5
    src.InsideOf = 1
    seti(src, "# Layout Rays", 1)
    seti(src, "# Analysis Rays", 1)
    setd(src, "Power(Watts)", 1.0)
    sd = src.SourcesData
    sd.RandomPolarization = False
    jx, jy = (1.0, 0.0) if pol == "TE" else (0.0, 1.0)
    sd.Jx, sd.Jy, sd.XPhase, sd.YPhase = jx, jy, 0.0, 0.0

    for i, zpos in ((4, -1.5), (5, 2.0)):
        det = make(i, T.DetectorRectangle)
        det.ZPosition = zpos
        det.Material = "ABSORB"
        setd(det, "X Half Width", 15.0)
        setd(det, "Y Half Width", 15.0)
        seti(det, "# X Pixels", 1)
        seti(det, "# Y Pixels", 1)
    return nce, holo, src


def trace(z, nce, split=False):
    rt = z.sys.Tools.OpenNSCRayTrace()
    rt.SplitNSCRays = bool(split)
    rt.ScatterNSCRays = False
    rt.UsePolarization = True
    rt.IgnoreErrors = True
    rt.SaveRays = False
    rt.ClearDetectors(0)
    rt.RunAndWaitForCompletion()
    rt.Close()
    out = []
    for i in (4, 5):
        ok, v = nce.GetDetectorData(i, 0, 0, 0.0)
        if not ok:
            raise RuntimeError(f"GetDetectorData failed on object {i}")
        out.append(v)
    return out


def main(name, pol="TE"):
    case = CASES[name]
    suf = "" if pol == "TE" else "_TM"
    z = Zos()
    print(f"OpticStudio {z.app.ZOSMajorVersion}.{z.app.ZOSMinorVersion}.{z.app.ZOSSPVersion}  licence {z.app.LicenseStatus}  pol {pol}", flush=True)
    rows = []
    failed = False
    try:
        nce, holo, src = build(z, case, pol=pol)
        os.makedirs(os.path.join(ROOT, "zemax"), exist_ok=True)
        zmx = os.path.join(ROOT, "zemax", f"kogelnik_native_{name}{suf}.zmx")
        z.sys.SaveAs(zmx)
        w1 = z.sys.SystemData.Wavelengths.GetWavelength(1)
        for tau in case["zemax_tilts_deg"]:
            src.TiltAboutX = float(tau)
            w1.Wavelength = case["lam_um"]
            pr, pt = trace(z, nce)
            rows.append(("angle", tau, case["lam_um"], pr, pt))
        src.TiltAboutX = float(case["zemax_tilt_bragg_deg"])
        for lam in case["zemax_lams_um"]:
            w1.Wavelength = float(lam)
            pr, pt = trace(z, nce)
            rows.append(("wavelength", case["zemax_tilt_bragg_deg"], lam, pr, pt))
        w1.Wavelength = case["lam_um"]
        z.sys.Save()
    except Exception:
        import traceback
        traceback.print_exc()    # os._exit below would otherwise swallow it
        failed = True
    finally:
        # A failed or empty run must never overwrite a good table.
        bad = failed or not rows
        out = os.path.join(ROOT, "results", f"zemax_{name}{suf}{'.FAILED' if bad else ''}.csv")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sweep", "source_tilt_x_deg", "lambda_um", "P_entry_side_W", "P_far_side_W"])
            w.writerows(rows)
        print(f"{len(rows)} rows -> {out}")
        try:
            z.close()
        except Exception as e:      # a dead server must not mask the real error
            print("close:", e)
        sys.stdout.flush()
        os._exit(1 if failed or not rows else 0)   # pythonnet + ZOS-API hangs on normal interpreter exit


if __name__ == "__main__":
    argv = sys.argv[1:]
    pol = "TE"
    if "--pol" in argv:
        i = argv.index("--pol")
        pol = argv[i + 1].upper()
        del argv[i:i + 2]
    main(argv[0], pol)
