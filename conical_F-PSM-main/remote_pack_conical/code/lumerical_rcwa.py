"""Drive the Ansys Lumerical RCWA solver on the article's own model
(``../lumerical/ansys_example/vhg_R_6um_model.fsp``, Knowledgebase article 37118679173267)
for the gratings in ``cases.py`` and write first-order efficiency tables in the pack's
common schema (the same schema the OpticStudio tables use, so ``compare.py`` overlays and
scores them with no change).

STATUS (2026-09-18): **VERIFIED at the solver** on Ansys Lumerical v261 (FDTD engine 8.35.4522),
Windows, Python 3.13.14.  All four benchmark cases ran and passed the gates of
EXECUTION_LUMERICAL.md section 6; every entry of the NAME MAP below was checked against this
build's own property and result lists (``results/lumerical_probe.txt``).  Everything that touches
the solver is still guarded: the driver refuses to configure anything until ``NAMES_VERIFIED``,
and writes ``*.FAILED.csv`` rather than overwriting a good table.  Order of work:

    1. python check_env.py                      -- Python, packages, lumapi found?
    2. python lumerical_rcwa.py --probe lumerical_R
                                                -- opens the model, DUMPS every property of
                                                   ::model, RCWA and the rectangle WITH ITS VALUE,
                                                   the model's setup script, and every available
                                                   result with its dataset structure
                                                   ->  results/lumerical_probe.txt
    3. check the NAME MAP below (``P``) against that printout; set ``NAMES_VERIFIED = True``
    4. python lumerical_rcwa.py --bragg <case>  -- one solve at Bragg: gate L0 (peak in the right
                                                   order), thickness read-back, energy
    5. python lumerical_rcwa.py --converge <case>
                                                -- harmonics and layers-per-period ladders at Bragg
                                                   and at the off-Bragg reference angles (gate L2)
    6. python lumerical_rcwa.py <case> --harmonics 3 --layers 160
                                                -- the two sweeps -> results/lumerical_rcwa_<case>.csv
    7. python compare.py                        -- overlays + scores

WHAT THE MODEL IS (read out of the .fsp, confirmed live 2026-09-18):
  * ``::model`` has user properties  "angle 1", "angle 2" [deg, IN THE MEDIUM], "holo type" ('R'|'T'),
    "n0", "Dn", "w_rec" [m];  its setup script builds  k1 = [sin a1, 0, cos a1],
    k2 = [sin a2, 0, +-cos a2] (minus for 'R'),  K = k1 - k2,  index = n0 + Dn cos((2 pi n0 / w_rec) K.r),
    x period dx = |w_rec / (n0 K_x)|, z period dz = |w_rec / (n0 K_z)|, sets the rectangle's "index"
    (a spatial equation string whose x, y, z are IN MICRONS), rectangle "z span" = dz,
    RCWA "x span" = dx, RCWA "z span" = 1.5 dz, RCWA background "index" = n0,
    RCWA "wavelength center" = w_rec.
  * That setup script RE-RUNS DURING EVERY SOLVE, so any geometry it owns that the driver
    overrides is silently restored unless the script is blanked first (see set_thickness) --
    ``solve`` therefore re-checks the thickness AFTER the run, not only before it.
  * The plate thickness is NOT a model property: it is the RCWA object's layer repetition
    (Interfaces tab), and "simulation z span" reports the whole region INCLUDING the background
    padding around the rectangle, so the grating itself is that minus (RCWA z span - rect z span).
  * There is no "layers" property: the layer boundaries are the RCWA object's interfaces, built
    from object faces by "interface reference positions" (see ``set_layers``).
  * There is no polarisation property either: s and p are solved together, and the TE quantity the
    article plots is result "grating_power", attribute Rs_power, order n = +1.
  * Lumerical's index equation is n = n0 + Dn cos(K.r), the note's definition.  Whether its RCWA
    keeps the Dn^2 terms of eps = n^2 was exactly what this leg had to tell; compare.py scores
    every Lumerical table against BOTH referee forms.  ANSWER (VERIFIED 2026-09-18): it keeps
    them -- rms against the exact referee is 1e-12 to 3e-5 on every sweep, against the
    eps = n0^2 + 2 n0 Dn cos form that the OpticStudio DLL implements it is 30-100x worse.
    See ../results/lumerical_verdict.md.

THE UNSLANTED TRAP.  ``unslanted_T`` has K_z = 0, so the model's dz is infinite and its setup
script cannot be used (the OpticStudio DLL crashes on this grating for the same reason).  Bypass,
implemented in ``configure_unslanted``: keep the model file but stop its setup script from
rebuilding the geometry (``P["model_script"]`` set to an empty string -- VERIFIED 2026-09-18 to
be accepted by lumapi, so no GUI step was needed), then set the rectangle's index equation to
n0 + Dn cos(2 pi/dx x) with dx = w_rec/(2 n0 sin a) IN MICRONS, the rectangle z span = the plate
thickness d, the RCWA x span = dx, z span >= d, and layer repetition "none".  If a future build
refuses to blank the script, do it ONCE in the GUI, save the file as
``../lumerical/unslanted_T_base.fsp`` and the driver will use that file instead.

CONVENTIONS.  cases.py: replay angle theta inside the medium, k = beta (sin theta, cos theta),
K = beta_r (k1 - k2).  The article's figure ("Rs power order +1 vs Theta") peaks at theta = 0 and
-28 deg for lumerical_R, which is cases.py's convention (theta_bragg = 0, second Bragg peak at
-28 deg).  Lumerical's own "theta" sign is NOT assumed: ``choose_sign`` compares the solved curve
with the referee under both signs and every candidate order, and refuses if neither fits.
VERIFIED 2026-09-18: the sign is +1, i.e. Lumerical's theta IS cases.py's theta, on all four cases.
Output column ``source_tilt_x_deg`` = -theta, as for the OpticStudio tables.

Usage
    python lumerical_rcwa.py --probe <case>      # writes results/lumerical_probe.txt, no table
    python lumerical_rcwa.py --bragg <case>      # gate L0, no table
    python lumerical_rcwa.py --converge <case>   # gate L2 ladders -> results/lumerical_converge_<case>.txt
    python lumerical_rcwa.py <case> [--harmonics N] [--layers N] [--full-structure]
    python lumerical_rcwa.py --bragg <case> --chi0 1e-3    # OPTIONAL absorption probe, one point (see below)

TM (added 2026-10-16, UNTESTED at the solver).  The solver computes s and p together, so every
run now ALSO tabulates the p (TM) attributes Rp_power / Tp_power -- the attribute names were
seen in the 2026-09-18 probe dump (results/lumerical_probe.txt, grating_power and total_energy
datasets) but never read -- and writes a second table in the common schema,
results/lumerical_rcwa_<case>_TM.csv, from the SAME solve (same order selection, same sign fit).
--probe now prints arrays in full (numpy's default print threshold truncated the 2026-09-18 dump
with ellipses).  --bragg and --converge print the TM values beside the TE ones, with the TM
referee.

ABSORPTION (optional, UNTESTED, --bragg only).  --chi0 X appends "-X i" to the rectangle's index
equation after the setup script has built it, i.e. n = (n0 - i X) + Dn cos(K.r).  Whether the
spatial index equation accepts a complex value is unknown; the gate is the energy sum: a solver
that ignores the imaginary part returns sum R+T = 1.000000 and the lossless value, and the driver
then raises.  Reference values: evidence/tm_reference_values.md (one absorbing point).
"""
import csv
import glob
import importlib.util
import math
import os
import sys

import numpy as np

import rcwa
from cases import CASES, deg

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FSP = os.path.join(ROOT, "lumerical", "ansys_example", "vhg_R_6um_model.fsp")
FSP_UNSLANTED = os.path.join(ROOT, "lumerical", "unslanted_T_base.fsp")     # optional, see the docstring
RES = os.path.join(ROOT, "results")

# Article-model inputs per case: (angle 1, angle 2, holo type), derived from cases.py:
# k1 = (sin a1, cos a1);  k2 = (sin a2, -cos a2) for 'R', (sin a2, +cos a2) for 'T'.
MODEL = {"lumerical_R": (0.0, 28.0, "R"), "lumerical_T": (0.0, 28.0, "T"),
         "zemax_kba_R": (-45.0, 0.0, "R"), "unslanted_T": (33.8627180249, -33.8627180249, "T"),
         "unslanted_T_rot5": (38.8627180249, -28.8627180249, "T")}

# Off-Bragg reference angles of the OpticStudio verdict (execution_original.md section 6) -- the
# discriminating points: a genuine RCWA lands on the RCWA column, not on Kogelnik's.
OFF_BRAGG = {"lumerical_R": (-14.0, 7.5), "lumerical_T": (5.0,), "zemax_kba_R": (-22.5,),
             "unslanted_T": (36.0,), "unslanted_T_rot5": (41.0,)}

# Settings held fixed while the other ladder runs (--converge, gate L2); generous on purpose.
H_LAYERS, L_HARMONICS = 80, 9
# Index-sampling cells along z per grating z-period.  100 cells over the shipped model's
# 1.5-z-period region (VERIFIED 2026-09-18 from the .fsp) = 100/1.5 per period.
MESH_Z_PER_PERIOD = 100.0 / 1.5

# ----------------------------------------------------------------------------------------------
# NAME MAP.  Every solver-facing name lives here.  Status of each entry:
#   VERIFIED 2026-09-18  the name appears in this build's own property/result list
#                        (Lumerical FDTD 8.35.4522, install v261) -- see results/lumerical_probe.txt
#   ASSUMED              still a guess
# Set NAMES_VERIFIED = True only after every entry has been checked against the probe printout.
# ----------------------------------------------------------------------------------------------
NAMES_VERIFIED = True
P = {
    # objects
    "model": "::model",                       # VERIFIED 2026-09-18
    "rcwa": "RCWA",                           # VERIFIED 2026-09-18
    "rect": "rectangle",                      # VERIFIED 2026-09-18
    # ::model user properties  (all VERIFIED 2026-09-18 in the ::model property list)
    "angle1": "angle 1",
    "angle2": "angle 2",
    "holotype": "holo type",                  # 'R' | 'T'
    "n0": "n0",
    "Dn": "Dn",                               # cosine AMPLITUDE: n = n0 + Dn cos
    "w_rec": "w_rec",                         # [m]
    "model_script": "setup script",           # VERIFIED 2026-09-18 (read back; blanking it is the unslanted bypass)
    # rectangle  (VERIFIED 2026-09-18)
    "rect_index": "index",                    # a spatial-equation string; ITS x,y,z ARE IN MICRONS
    "rect_zspan": "z span",
    # RCWA geometry / background  (VERIFIED 2026-09-18)
    "rcwa_xspan": "x span",
    "rcwa_zspan": "z span",
    "rcwa_index": "index",                    # background index (1.5)
    # RCWA layer repetition (article 36965778575379)  (VERIFIED 2026-09-18)
    "rep_mode": "layer repetition",           # "none" | "multiple" | "thickness"
    "rep_start": "layer repetition start interface",   # indexes a ROW of "interface reference positions"
    "rep_end": "layer repetition end interface",
    "rep_thickness": "layer repetition thickness",     # metres
    "rep_multiple": "layer repetition multiple",       # read back: 30.8119 periods for lumerical_R
    "sim_zspan": "simulation z span",         # read-only; INCLUDES the RCWA region's background padding
    # RCWA slicing into homogeneous layers.  There is no "layers" property: the layer boundaries are
    # the RCWA object's INTERFACES, and the shipped file generates them from the rectangle's faces:
    #   interface reference positions = [["::model::rectangle","max",10],["::model::rectangle","min",1]]
    # VERIFIED 2026-09-18 by experiment: the "max" row's third column N yields N+1 evenly spaced
    # "interface actual positions" across the rectangle's z span, i.e. N layers per z-period
    # (N = 4 -> 5 positions, N = 20 -> 21 positions, span always the rectangle's z span).
    "interfaces_ref": "interface reference positions",    # VERIFIED 2026-09-18
    "interfaces_actual": "interface actual positions",    # VERIFIED 2026-09-18 (read-only)
    # RCWA harmonics (article 42661761082899)  (VERIFIED 2026-09-18)
    "k_domain": "k vectors domain",           # "circular" (shipped) | "rectangular"
    "k_max": "max number k vectors",          # max |order| when the domain is circular; 50 shipped
    "ku_max": "max number ku",
    "kv_max": "max number kv",
    # RCWA excitation.  VERIFIED 2026-09-18: this object has NO polarisation property -- it solves
    # s and p together and reports them as separate result attributes (Rs_power / Rp_power).
    # TE ("S") is therefore just the Rs_/Ts_ attributes; nothing to set.
    "theta_mode": "incident angle",           # VERIFIED 2026-09-18: "range" (shipped) | "table"
    "theta_table": "incident angle table",    # VERIFIED 2026-09-18 (N x 2, columns theta and phi, degrees)
    "theta_actual": "incident angles actual", # VERIFIED 2026-09-18 (read-only N x 2: what will be solved)
    "theta_min": "minimum theta",             # VERIFIED 2026-09-18
    "theta_max": "maximum theta",             # VERIFIED 2026-09-18
    "theta_n": "theta points",                # VERIFIED 2026-09-18
    "theta_spacing": "theta spacing",         # VERIFIED 2026-09-18 ("linear")
    "phi_min": "minimum phi",                 # VERIFIED 2026-09-18 (0: plane of incidence is x-z)
    "phi_max": "maximum phi",                 # VERIFIED 2026-09-18
    "phi_n": "phi points",                    # VERIFIED 2026-09-18
    "wl_center": "wavelength center",         # VERIFIED 2026-09-18
    "wl_span": "wavelength span",             # VERIFIED 2026-09-18
    "use_wl_spacing": "use wavelength spacing",   # VERIFIED 2026-09-18 (1: sample linearly in lambda)
    "freq_points": "frequency points",        # VERIFIED 2026-09-18
    "prop_axis": "propagation axis",          # VERIFIED 2026-09-18 ("z")
    "prop_dir": "propagation direction",      # VERIFIED 2026-09-18 ("both": forward AND backward incidence)
    "mesh_z": "mesh cells z",                 # VERIFIED 2026-09-18 (index sampling along z over the WHOLE region)
    "mesh_z_actual": "actual mesh cells z",   # VERIFIED 2026-09-18 (read-only)
    # results.  VERIFIED 2026-09-18: the available results are x,y,z,num_k,peak_memory,
    # simulation_run_time,total_threads,total_energy,substrate,grating_orders,grating_power,
    # grating_characterization.  "grating_characterization" carries complex AMPLITUDES (Rss...),
    # "grating_orders" the efficiencies of every retained order including evanescent ones, and
    # "grating_power" the POWER in each PROPAGATING order -- the article's figure is
    # "Rs power order +1 vs Theta", so grating_power/Rs_power is the quantity to reproduce.
    # Dataset layout: [lambda, theta, n, m, forward|backward]; m is the y order (single 0 here).
    "result_grating": "grating_power",        # VERIFIED 2026-09-18
    "res_R": "Rs_power",                      # VERIFIED 2026-09-18 (s = TE, reflected)
    "res_T": "Ts_power",                      # VERIFIED 2026-09-18
    "res_Rp": "Rp_power",                     # name VERIFIED 2026-09-18 (grating_power attribute list in the probe dump); values UNTESTED (never read)
    "res_Tp": "Tp_power",                     # idem
    "res_orders": "n",                        # VERIFIED 2026-09-18 (the x diffraction order)
    "res_m": "m",                             # VERIFIED 2026-09-18 (the y order; size 1)
    "res_theta": "theta",                     # VERIFIED 2026-09-18
    "res_lambda": "lambda",                   # VERIFIED 2026-09-18 ("f" is handled as a fallback)
    "res_fwd": "forward",                     # VERIFIED 2026-09-18 (selector [1,0] on the last axis)
    "result_energy": "total_energy",          # VERIFIED 2026-09-18 (a SEPARATE result)
    "res_energy_R": "Rs",                     # VERIFIED 2026-09-18
    "res_energy_T": "Ts",                     # VERIFIED 2026-09-18
    "res_energy_Rp": "Rp",                    # name VERIFIED 2026-09-18 (total_energy keys Rs, Ts, Rp, Tp in the probe dump); values UNTESTED
    "res_energy_Tp": "Tp",                    # idem
}


# ----------------------------------------------------------------------------- lumapi ---------
def load_lumapi():
    """Resolve lumapi the way the OpticStudio driver resolves ZOSAPI: env var, then a scan.  Never a literal."""
    cands = [os.environ.get("LUMAPI_PATH", "")]
    for pat in (r"C:\Program Files\Lumerical\v*\api\python\lumapi.py",
                r"C:\Program Files\ANSYS Inc\v*\Lumerical\api\python\lumapi.py",
                r"C:\Program Files\Ansys Inc\v*\Optical Products\Lumerical\api\python\lumapi.py",
                r"C:\Program Files\ANSYS Inc\v*\Optics\Lumerical\api\python\lumapi.py"):
        cands += sorted(glob.glob(pat), reverse=True)
    for c in cands:
        if c and os.path.isfile(c):
            api_dir = os.path.dirname(c)
            root = os.path.dirname(os.path.dirname(api_dir))
            for d in (os.path.join(root, "bin"), api_dir):
                if os.path.isdir(d):
                    try:
                        os.add_dll_directory(d)
                    except Exception:
                        pass
            sys.path.insert(0, api_dir)
            spec = importlib.util.spec_from_file_location("lumapi", c)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            print("lumapi:", c, flush=True)
            return mod
    raise RuntimeError("lumapi.py not found -- set LUMAPI_PATH to its full path")


def open_model(lumapi, path=FSP):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    return lumapi.FDTD(filename=path, hide=True)


# ----------------------------------------------------------------------------- probe ----------
def try_(label, fn, out):
    try:
        v = fn()
        out.append(f"--- {label}\n{v}\n")
        return v
    except Exception as e:
        out.append(f"--- {label}: FAILED: {str(e).splitlines()[0] if str(e) else type(e).__name__}\n")
        return None


def describe(tag, obj, out):
    """Summarise a result dataset: every key, its shape, dtype and (if small) its values."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            a = np.asarray(v)
            out.append(f"    {tag} key {k!r}: shape {a.shape} dtype {a.dtype}"
                       + (f" values {a.ravel()[:12]}" if a.size <= 64 else "") + "\n")
    else:
        a = np.asarray(obj)
        out.append(f"    {tag}: shape {a.shape} dtype {a.dtype}"
                   + (f" values {a.ravel()[:12]}" if a.size <= 64 else "") + "\n")


def probe(fdtd, name):
    """Dump what this Lumerical version really exposes.  Nothing here is ASSUMED to work: each
    call is tried independently and its failure is recorded, so the printout is always complete.

    Three passes per object: the property list, every property in it read back with its real
    value (this is what settles the ASSUMED names in ``P``), and then every property ``P`` names.
    Finally every available result is fetched and its dataset structure described."""
    out = [f"Lumerical probe, case {name}, file {FSP}\n"]
    # Print every array in full: the 2026-09-18 dump was elided by numpy's default threshold, which
    # left the Rp_power / Tp_power arrays as fragments (2026-10-16).
    np.set_printoptions(threshold=sys.maxsize, linewidth=160)
    try_("version", lambda: fdtd.version(), out)
    for obj in (P["model"], P["rcwa"], P["rect"]):
        # 1. getnamed with the object name only: in Lumerical script this returns the property list
        lst = try_(f"getnamed({obj!r}) property list", lambda o=obj: fdtd.getnamed(o), out)
        # 2. every property that list names, read back with its value
        out.append(f"--- {obj!r}: every property with its current value\n")
        for prop in sorted(set(str(lst or "").split("\n"))):
            prop = prop.strip()
            if not prop:
                continue
            try:
                v = fdtd.getnamed(obj, prop)
                a = np.asarray(v)
                s = repr(v) if a.size <= 12 else f"<shape {a.shape} dtype {a.dtype}>"
                out.append(f"    {prop!r} = {s[:400]}\n")
            except Exception as e:
                out.append(f"    {prop!r}: unreadable ({str(e).splitlines()[0][:120]})\n")
        # 3. every property we intend to touch, read back
        for key, prop in P.items():
            if key in ("model", "rcwa", "rect") or key.startswith("res_") or key.startswith("result_"):
                continue
            try_(f"getnamed({obj!r}, {prop!r})   [P[{key!r}]]", lambda o=obj, p=prop: fdtd.getnamed(o, p), out)
    try_("model setup script", lambda: fdtd.getnamed(P["model"], P["model_script"]), out)
    names = try_("available results of RCWA (getresult with no result name)",
                 lambda: fdtd.getresult(P["rcwa"]), out)
    for rn in sorted(set(str(names or "").split("\n"))):
        rn = rn.strip()
        if not rn:
            continue
        r = try_(f"getresult(RCWA, {rn!r})", lambda n=rn: fdtd.getresult(P["rcwa"], n), out)
        if r is not None:
            describe(rn, r, out)
    txt = "".join(out)
    path = os.path.join(RES, "lumerical_probe.txt")
    os.makedirs(RES, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)
    print(txt)
    print(f"-> {path}", flush=True)


# ----------------------------------------------------------------------------- configure ------
def require_names():
    if not NAMES_VERIFIED:
        raise RuntimeError("NAME MAP not verified: run --probe, fix every ASSUMED entry in P, set NAMES_VERIFIED = True")


_SHIPPED_SCRIPT = {}


CHI0 = 0.0      # optional absorption probe (--chi0), applied by set_model; 0 = the benchmark's lossless grating


def set_model(fdtd, case, name):
    a1, a2, typ = MODEL[name]
    fdtd.switchtolayout()
    if "s" in _SHIPPED_SCRIPT:
        fdtd.setnamed(P["model"], P["model_script"], _SHIPPED_SCRIPT["s"])   # --full-structure blanks it
    else:
        _SHIPPED_SCRIPT["s"] = fdtd.getnamed(P["model"], P["model_script"])
    for key, val in (("angle1", a1), ("angle2", a2), ("holotype", typ), ("n0", case["n0"]),
                     ("Dn", case["n1"]), ("w_rec", case["lam_r_um"] * 1e-6)):
        fdtd.setnamed(P["model"], P[key], val)
    fdtd.runsetup()          # force ::model's setup script to rebuild the rectangle and the RCWA region
    if CHI0 > 0.0:
        # UNTESTED 2026-10-16: complex index through the rectangle's spatial index equation.  The setup
        # script re-runs at every solve and would restore the real equation, so it is blanked here as
        # for --full-structure; set_model restores it on the next call.
        expr = str(fdtd.getnamed(P["rect"], P["rect_index"]))
        fdtd.setnamed(P["model"], P["model_script"], "")
        fdtd.setnamed(P["rect"], P["rect_index"], f"({expr})-{CHI0:g}i")
        print(f"   absorption probe: rectangle index = ({expr})-{CHI0:g}i   (UNTESTED: gate = energy sum < 1)", flush=True)


def set_thickness(fdtd, case, dz, full_structure=False):
    """Plate thickness d through layer repetition (mode 'thickness'), or by building the full
    structure (rectangle z span = d, repetition off) when asked.

    VERIFIED 2026-09-18: ``mesh cells z`` samples the index over the WHOLE RCWA region, so the
    shipped value (100 cells over the 1.5 z-periods the repetition case simulates) would alias the
    sinusoid once the region becomes the whole plate.  It is therefore always set from the
    region's size in z-periods, which reproduces the shipped value in the repetition case."""
    if full_structure:
        # The ::model setup script owns "rectangle z span" and re-runs on every solve, so it would
        # undo the resize between configure() and the result.  VERIFIED 2026-09-18: left in place
        # it does exactly that, and the full structure returns 0.002236 at Bragg -- which is
        # tanh^2 of ONE z-period (0.1947 um), not of 6 um.  It has already done its job (the index
        # equation and dz) by the time we get here, so blank it; set_model restores it.
        fdtd.setnamed(P["model"], P["model_script"], "")
        fdtd.setnamed(P["rcwa"], P["rep_mode"], "none")
        fdtd.setnamed(P["rect"], P["rect_zspan"], case["d_um"] * 1e-6)
        fdtd.setnamed(P["rcwa"], P["rcwa_zspan"], case["d_um"] * 1e-6 * 1.5)
    else:
        fdtd.setnamed(P["rcwa"], P["rep_mode"], "thickness")
        # The shipped file already repeats the rectangle's interfaces; the start/end interface
        # indices are whatever the file holds -- read them back rather than guess.
        print("   layer repetition interfaces:", fdtd.getnamed(P["rcwa"], P["rep_start"]), "->",
              fdtd.getnamed(P["rcwa"], P["rep_end"]), flush=True)
        fdtd.setnamed(P["rcwa"], P["rep_thickness"], np.array([[case["d_um"] * 1e-6]]))
    cells = MESH_Z_PER_PERIOD * float(fdtd.getnamed(P["rcwa"], P["rcwa_zspan"])) / dz
    fdtd.setnamed(P["rcwa"], P["mesh_z"], int(math.ceil(cells)))


def read_thickness(fdtd, full_structure=False):
    """The grating thickness the solver will really simulate, in metres (gate L0/L2: it must equal d,
    not a whole-period truncation).

    VERIFIED 2026-09-18: ``simulation z span`` is the span of the whole RCWA region after the
    repetition has stretched it, so it includes the background padding the shipped model leaves
    around the rectangle (RCWA z span = 1.5 dz against a rectangle z span of dz).  The grating
    itself is therefore ``simulation z span - (RCWA z span - rectangle z span)``: for the shipped
    file 6.097365074858495e-6 - (2.9209522457548653e-7 - 1.94730149716991e-7) = 6.000000e-6 m.
    The padding is background index n0, which is exactly the benchmark's no-Fresnel-face condition."""
    if full_structure or str(fdtd.getnamed(P["rcwa"], P["rep_mode"])) == "none":
        return float(fdtd.getnamed(P["rect"], P["rect_zspan"]))
    pad = float(fdtd.getnamed(P["rcwa"], P["rcwa_zspan"])) - float(fdtd.getnamed(P["rect"], P["rect_zspan"]))
    return float(fdtd.getnamed(P["rcwa"], P["sim_zspan"])) - pad


def set_excitation(fdtd, thetas_deg, lams_m):
    """The RCWA object sweeps theta either as min/max/points ("range") or as an explicit
    (theta, phi) table.  The table is used whenever the requested angles are not a uniform
    grid -- the convergence ladders ask for three specific angles, and cases.py's sweep grids
    must be reproduced exactly so compare.py needs no interpolation.

    There is no polarisation property: s and p are solved together and reported as separate
    result attributes (VERIFIED 2026-09-18)."""
    th, lam = np.asarray(thetas_deg, float), np.asarray(lams_m, float)
    uniform = th.size < 3 or np.allclose(np.diff(th), th[1] - th[0])
    if uniform and th.size > 1:
        fdtd.setnamed(P["rcwa"], P["theta_mode"], "range")
        fdtd.setnamed(P["rcwa"], P["theta_spacing"], "linear")
        # "minimum theta" is rejected if it would exceed the CURRENT "maximum theta" (VERIFIED
        # 2026-09-18: unslanted_T asks for 33.86 deg while the shipped file still holds 30).
        # Widen the window first, then narrow it, so the two never cross.
        lo = min(float(th.min()), float(fdtd.getnamed(P["rcwa"], P["theta_min"])))
        hi = max(float(th.max()), float(fdtd.getnamed(P["rcwa"], P["theta_max"])))
        fdtd.setnamed(P["rcwa"], P["theta_min"], lo)
        fdtd.setnamed(P["rcwa"], P["theta_max"], hi)
        fdtd.setnamed(P["rcwa"], P["theta_min"], float(th.min()))
        fdtd.setnamed(P["rcwa"], P["theta_max"], float(th.max()))
        fdtd.setnamed(P["rcwa"], P["theta_n"], int(th.size))
        fdtd.setnamed(P["rcwa"], P["phi_min"], 0.0)
        fdtd.setnamed(P["rcwa"], P["phi_max"], 0.0)
        fdtd.setnamed(P["rcwa"], P["phi_n"], 1)
    else:
        fdtd.setnamed(P["rcwa"], P["theta_mode"], "table")
        fdtd.setnamed(P["rcwa"], P["theta_table"], np.column_stack([th, np.zeros_like(th)]))
    got = np.asarray(fdtd.getnamed(P["rcwa"], P["theta_actual"]))
    if got.shape[0] != th.size or not np.allclose(np.sort(got[:, 0]), np.sort(th), atol=1e-9):
        raise RuntimeError(f"asked {th.size} angles, solver will use {got.shape[0]}: {got[:, 0][:6]} ...")
    fdtd.setnamed(P["rcwa"], P["use_wl_spacing"], 1)
    fdtd.setnamed(P["rcwa"], P["wl_center"], float(0.5 * (lam.min() + lam.max())))
    fdtd.setnamed(P["rcwa"], P["wl_span"], float(lam.max() - lam.min()))
    fdtd.setnamed(P["rcwa"], P["freq_points"], int(lam.size))


def set_layers(fdtd, layers, periods=1.0):
    """Slice the graded-index rectangle into homogeneous layers.

    VERIFIED 2026-09-18: the RCWA object has no "layers" property.  Its layer boundaries are the
    INTERFACES, generated from object faces by ``interface reference positions``, a table of
    [object, face, count] rows.  The shipped file holds
    [["::model::rectangle","max",10],["::model::rectangle","min",1]] and reports 11 evenly spaced
    ``interface actual positions`` across the rectangle's z span -- 10 layers per z-period, the
    number the article quotes.  Setting the "max" row's count to N gives N+1 positions (checked at
    N = 4 and N = 20), and ``layer repetition start/end interface`` = 1/2 index the ROWS of this
    table, so the repetition still spans the whole rectangle whatever N is.

    ``periods`` > 1 is the --full-structure case, where the rectangle is the whole plate rather
    than one z-period, so the same layers-per-period needs that many times more interfaces."""
    n = int(round(layers * periods))
    fdtd.setnamed(P["rcwa"], P["interfaces_ref"],
                  [[f"::model::{P['rect']}", "max", float(n)], [f"::model::{P['rect']}", "min", 1.0]])
    got = np.asarray(fdtd.getnamed(P["rcwa"], P["interfaces_actual"])).size
    if got != n + 1:
        raise RuntimeError(f"asked {n} layers, solver made {got} interfaces (expected {n + 1})")


def set_numerics(fdtd, harmonics=None, layers=None, periods=1.0):
    if harmonics is not None:
        # 1-D grating along x, simulation region "2D Y-normal": a rectangular k domain with
        # ku = harmonics retains the orders -harmonics .. +harmonics along x and none along y.
        fdtd.setnamed(P["rcwa"], P["k_domain"], "rectangular")
        fdtd.setnamed(P["rcwa"], P["ku_max"], int(harmonics))
        fdtd.setnamed(P["rcwa"], P["kv_max"], 0)
    if layers is not None:
        set_layers(fdtd, layers, periods)


def configure(fdtd, case, name, thetas_deg, lams_m, harmonics=None, layers=None, full_structure=False):
    require_names()
    if name == "unslanted_T":
        return configure_unslanted(fdtd, case, thetas_deg, lams_m, harmonics, layers)
    set_model(fdtd, case, name)
    dz = float(fdtd.getnamed(P["rect"], P["rect_zspan"]))          # one z-period, from the setup script
    set_thickness(fdtd, case, dz, full_structure)
    set_excitation(fdtd, thetas_deg, lams_m)
    set_numerics(fdtd, harmonics, layers, case["d_um"] * 1e-6 / dz if full_structure else 1.0)
    d_sim = read_thickness(fdtd, full_structure)
    print(f"   z-period {dz * 1e6:.6f} um; simulated thickness {d_sim * 1e6:.4f} um (asked {case['d_um']:.4f}); "
          f"mesh cells z {fdtd.getnamed(P['rcwa'], P['mesh_z'])!r}", flush=True)
    if abs(d_sim - case["d_um"] * 1e-6) > 5e-9:
        raise RuntimeError(f"thickness mismatch: solver will simulate {d_sim * 1e6:.4f} um, asked {case['d_um']} um "
                           "(whole-period truncation? use --full-structure and compare)")


def configure_unslanted(fdtd, case, thetas_deg, lams_m, harmonics=None, layers=None):
    """K_z = 0 bypass (see the module docstring).  Everything here is the same grating cases.py
    defines: n = n0 + n1 cos(2 pi x / dx), dx = lam_r / (2 n0 sin a), plate thickness d.

    VERIFIED 2026-09-18: the rectangle's index equation takes x, y, z in MICRONS, not metres --
    the shipped setup script writes ``(2*pi*1.5/0.55)*(x*(-0.469472)+...)`` with w_rec in um."""
    a = MODEL["unslanted_T"][0] * deg
    dx_um = case["lam_r_um"] / (2.0 * case["n0"] * math.sin(a))
    fdtd.switchtolayout()
    if not os.path.isfile(FSP_UNSLANTED):
        fdtd.setnamed(P["model"], P["model_script"], "")             # stop the script rebuilding dz = inf
    fdtd.setnamed(P["rect"], P["rect_index"], f"{case['n0']}+{case['n1']}*cos(2*pi*x/{dx_um})")
    fdtd.setnamed(P["rect"], P["rect_zspan"], case["d_um"] * 1e-6)
    fdtd.setnamed(P["rcwa"], P["rcwa_xspan"], dx_um * 1e-6)
    fdtd.setnamed(P["rcwa"], P["rcwa_zspan"], case["d_um"] * 1e-6 * 1.5)
    fdtd.setnamed(P["rcwa"], P["rcwa_index"], case["n0"])
    fdtd.setnamed(P["rcwa"], P["rep_mode"], "none")
    # n(x) only: nothing to resolve along z, so the shipped z sampling is ample (it must still be
    # set, because a previous --full-structure configure in the same session raises it).
    fdtd.setnamed(P["rcwa"], P["mesh_z"], 100)
    set_excitation(fdtd, thetas_deg, lams_m)
    set_numerics(fdtd, harmonics, layers)
    print(f"   unslanted bypass: dx = {dx_um:.6f} um, rectangle z span = {case['d_um']} um, repetition off", flush=True)


# ----------------------------------------------------------------------------- solve ----------
def run(fdtd):
    try:
        fdtd.run(P["rcwa"])          # VERIFIED 2026-09-18: this build accepts a named solver
    except Exception:
        fdtd.run()                   # fallback: run the active solver


def solve(fdtd, case, full_structure=False):
    """Run and read the grating powers, re-checking the geometry AFTERWARDS: the solve itself
    re-executes ::model's setup script, so a thickness that was right at configure time can still
    be wrong at the result (VERIFIED 2026-09-18, --full-structure)."""
    run(fdtd)
    d_sim = read_thickness(fdtd, full_structure)
    if abs(d_sim - case["d_um"] * 1e-6) > 5e-9:
        raise RuntimeError(f"the run reset the geometry: it simulated {d_sim * 1e6:.4f} um, "
                           f"not the {case['d_um']} um asked")
    res = fdtd.getresult(P["rcwa"], P["result_grating"])
    if not isinstance(res, dict):
        raise RuntimeError(f"{P['result_grating']} is not a dict; adapt orders_array()")
    return res


def axis(res, key, fallback=None):
    if key in res:
        return np.real(np.asarray(res[key])).ravel()
    if fallback and fallback in res:
        return 299792458.0 / np.real(np.asarray(res[fallback])).ravel()
    raise KeyError(f"result has no {key!r} (keys: {list(res)})")


def orders_array(res, reflection, pol="TE"):
    """[order, freq, theta] real efficiencies for FORWARD incidence; TE reads the s attributes
    (Rs_power / Ts_power), TM the p attributes (Rp_power / Tp_power; names seen in the probe dump,
    values UNTESTED 2026-10-16).

    VERIFIED 2026-09-18: ``grating_power``'s attributes are shaped
    [lambda, theta, n, m, forward|backward] -- n the x diffraction order, m the y order (a single
    0 for this 2D Y-normal region), and the last axis selected by the dataset's own ``forward`` =
    [1, 0] / ``backward`` = [0, 1] vectors, i.e. which face the plane wave enters by."""
    kR, kT = ("res_R", "res_T") if pol == "TE" else ("res_Rp", "res_Tp")
    arr = np.real(np.asarray(res[P[kR] if reflection else P[kT]]))
    n_ord = np.asarray(res[P["res_orders"]]).size
    n_th = axis(res, P["res_theta"]).size
    n_f = axis(res, P["res_lambda"], fallback="f").size
    n_m = np.asarray(res[P["res_m"]]).size
    if n_m != 1:
        raise RuntimeError(f"the y order axis has {n_m} entries; this benchmark is 1-D in x")
    if arr.shape != (n_f, n_th, n_ord, n_m, 2):
        raise RuntimeError(f"unexpected {P[kR] if reflection else P[kT]} shape {arr.shape}; "
                           f"expected (lambda {n_f}, theta {n_th}, order {n_ord}, m 1, 2)")
    fwd = int(np.argmax(np.asarray(res[P["res_fwd"]]).ravel()))
    return np.transpose(arr[:, :, :, 0, fwd], (2, 0, 1))


def total_energy(fdtd, reflection, pol="TE"):
    """Sum over every propagating order of the entry-side / far-side power, forward incidence;
    TE from Rs/Ts, TM from Rp/Tp (TM names seen in the probe dump, values UNTESTED 2026-10-16)."""
    te = fdtd.getresult(P["rcwa"], P["result_energy"])
    fwd = int(np.argmax(np.asarray(te[P["res_fwd"]]).ravel()))
    kR, kT = ("res_energy_R", "res_energy_T") if pol == "TE" else ("res_energy_Rp", "res_energy_Tp")
    r = np.real(np.asarray(te[P[kR]]))[..., fwd]
    t = np.real(np.asarray(te[P[kT]]))[..., fwd]
    return (r, t) if reflection else (t, r)


def referee(case, thetas_deg, lam_nm, reflection, eps_form="exact", pol="TE", chi0=0.0):
    k0 = 2 * np.pi / np.asarray(lam_nm, float)
    th = np.asarray(thetas_deg, float)
    psi = np.where(np.cos(th * deg + case["psi"]) >= 0, case["psi"], case["psi"] + np.pi)   # compare.py's K.k rule
    eR, eT = rcwa.signal_efficiency(k0, case["n0"], case["n1"], case["K"], psi, th * deg + psi, case["d_nm"],
                                    L=10, eps_form=eps_form, pol=pol, chi0=chi0)
    return eR if reflection else eT


def conj_index(orders, k):
    """Index of the conjugate integer order -n of order n, or None if the solver did not keep it."""
    hit = np.flatnonzero(orders == -orders[k])
    return int(hit[0]) if hit.size else None


def first_order_curve(eta_by_order, orders, k, thetas_deg, case):
    """The "first order" as compare.py labels it: the Bragg integer order while cos(theta + psi) >= 0,
    its conjugate beyond.  For a strongly slanted transmission grating (lumerical_T, psi = 76 deg)
    the swap happens inside the sweep, at theta = 14 deg."""
    kc = conj_index(orders, k)
    if kc is None:
        return eta_by_order[k]
    flip = np.cos(np.asarray(thetas_deg, float) * deg + case["psi"]) < 0
    use = np.where(flip, kc, k)
    return eta_by_order[use, np.arange(eta_by_order.shape[1])]


def choose_sign(theta_lum, eta_by_order, orders, case, lam_nm, reflection):
    """Lumerical's theta sign is not assumed.  For each sign and each candidate Bragg order, the rms
    distance to the referee; pick the best; refuse unless it is a genuine fit.

    The candidate curve is built with the SAME K.k order-swap rule the referee and the output rows
    use.  VERIFIED 2026-09-18: comparing one fixed integer order across the whole sweep instead
    makes the search ill-posed wherever the swap falls inside the sweep -- on lumerical_T it chose
    the wrong theta sign and reported rms 2.6e-2, against 2.3e-5 for the cases with no swap."""
    best = None
    for sign in (+1.0, -1.0):
        th = sign * theta_lum
        ref = referee(case, th, lam_nm, reflection)
        m = np.isfinite(ref)
        for k in range(eta_by_order.shape[0]):
            e = first_order_curve(eta_by_order, orders, k, th, case)
            rms = float(np.sqrt(np.mean((e[m] - ref[m]) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, sign, k)
    rms, sign, k = best
    print(f"   theta sign {sign:+.0f}, Bragg order {orders[k]:+.0f} (index {k}): rms vs referee {rms:.2e}", flush=True)
    if rms > 0.05:
        raise RuntimeError(f"no (sign, order) combination fits the referee (best rms {rms:.2e}): geometry or result layout is wrong")
    return sign, k


# ----------------------------------------------------------------------------- modes ----------
def bragg_check(fdtd, case, name, harmonics, layers, full):
    reflection = MODEL[name][2] == "R"
    thb = case["theta_bragg_deg"]
    configure(fdtd, case, name, [thb], [case["lam_um"] * 1e-6], harmonics, layers, full)
    res = solve(fdtd, case, full)
    arr = orders_array(res, reflection)[:, 0, 0]
    orders = np.asarray(res[P["res_orders"]]).ravel()
    k = int(np.argmax(arr))
    other = orders_array(res, not reflection)[:, 0, 0]
    print(f"   orders {orders.tolist()}\n   {'R' if reflection else 'T'} by order: {np.round(arr, 6).tolist()}\n"
          f"   {'T' if reflection else 'R'} by order: {np.round(other, 6).tolist()}", flush=True)
    tot = float(arr.sum() + other.sum())
    er, et = total_energy(fdtd, reflection)
    tot_e = float(np.asarray(er).ravel()[0] + np.asarray(et).ravel()[0])
    print(f"   Bragg peak in order index {k} (order {orders[k]}): eta = {arr[k]:.6f};  sum over orders R+T = {tot:.6f}"
          f";  total_energy R+T = {tot_e:.6f}", flush=True)
    ref = float(referee(case, [thb], case["lam_nm"], reflection, chi0=CHI0)[0])
    refl = float(referee(case, [thb], case["lam_nm"], reflection, "linear", chi0=CHI0)[0])
    print(f"   referee at Bragg (TE): exact {ref:.6f}   eps-linear {refl:.6f}   (execution_original.md section 6)", flush=True)
    # TM from the same solve (2026-10-16; UNTESTED at the solver)
    arr_p = orders_array(res, reflection, "TM")[:, 0, 0]
    other_p = orders_array(res, not reflection, "TM")[:, 0, 0]
    er_p, et_p = total_energy(fdtd, reflection, "TM")
    tot_p = float(np.asarray(er_p).ravel()[0] + np.asarray(et_p).ravel()[0])
    ref_tm = float(referee(case, [thb], case["lam_nm"], reflection, pol="TM", chi0=CHI0)[0])
    refl_tm = float(referee(case, [thb], case["lam_nm"], reflection, "linear", pol="TM", chi0=CHI0)[0])
    print(f"   TM (p) by order: {np.round(arr_p, 6).tolist()}; peak in order index {int(np.argmax(arr_p))}: {arr_p.max():.6f};"
          f" total_energy R+T = {tot_p:.6f}", flush=True)
    print(f"   referee at Bragg (TM): exact {ref_tm:.6f}   eps-linear {refl_tm:.6f}   (evidence/tm_reference_values.md)", flush=True)
    if CHI0 > 0.0:
        if tot_e > 1.0 - 1e-4:
            raise RuntimeError(f"absorption probe: energy sum {tot_e:.6f} is still 1 -- the imaginary index was ignored")
        print(f"   absorption probe: TE sum {tot_e:.6f} / TM sum {tot_p:.6f}; referee sums for chi0={CHI0:g}: see evidence/tm_reference_values.md", flush=True)
    elif abs(tot_e - 1.0) > 1e-3:
        raise RuntimeError(f"energy not conserved: R + T = {tot_e:.6f} (gate L0/Z3)")
    if arr[k] < 0.3:
        raise RuntimeError(f"no order carries the Bragg peak ({arr[k]:.4f}); geometry mapping is wrong")
    return k


def converge(fdtd, case, name, full):
    """Gate L2: a ladder in harmonics at a fixed generous layer count, and a ladder in layers per
    z-period at a fixed generous harmonic count.  Both ladders set BOTH numbers explicitly, so no
    rung inherits the previous one's settings."""
    reflection = MODEL[name][2] == "R"
    thetas = [case["theta_bragg_deg"]] + list(OFF_BRAGG[name])
    lines = [f"convergence, {name}, thetas {thetas} deg, lambda {case['lam_um']} um\n"]
    ref = referee(case, thetas, case["lam_nm"], reflection)
    refl = referee(case, thetas, case["lam_nm"], reflection, "linear")
    ref_tm = referee(case, thetas, case["lam_nm"], reflection, pol="TM")
    refl_tm = referee(case, thetas, case["lam_nm"], reflection, "linear", pol="TM")
    lines.append(f"referee exact   {np.round(ref, 6).tolist()}\nreferee linear  {np.round(refl, 6).tolist()}\n"
                 f"referee exact  TM {np.round(ref_tm, 6).tolist()}\nreferee linear TM {np.round(refl_tm, 6).tolist()}\n")
    for label, rungs in (("harmonics", [(h, H_LAYERS) for h in (3, 5, 9, 13)]),
                         ("layers per period", [(L_HARMONICS, l) for l in (10, 20, 40, 80, 160)])):
        for h, l in rungs:
            try:
                configure(fdtd, case, name, thetas, [case["lam_um"] * 1e-6], h, l, full)
                res = solve(fdtd, case, full)
                arr = orders_array(res, reflection)[:, 0, :]
                arr_p = orders_array(res, reflection, "TM")[:, 0, :]
                k = int(np.argmax(arr[:, 0]))
                th = axis(res, P["res_theta"])
                o = np.argsort(th)
                vals, vals_p = arr[k], arr_p[k]
                if not np.allclose(th, thetas):
                    vals = np.interp(thetas, th[o], vals[o])
                    vals_p = np.interp(thetas, th[o], vals_p[o])
                lines.append(f"{label} h={h} l={l}: {np.round(vals, 6).tolist()}   TM: {np.round(vals_p, 6).tolist()}\n")
            except Exception as e:
                lines.append(f"{label} h={h} l={l}: FAILED {str(e).splitlines()[0]}\n")
            print(lines[-1], end="", flush=True)
    path = os.path.join(RES, f"lumerical_converge_{name}.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"-> {path}  (gate L2: consecutive rungs must differ by < 1e-4)", flush=True)


def _rows_from(arr, other, reflection, orders, k, kc, sign, th_l, case, harmonics, layers):
    """Angle-sweep rows in the common schema for one polarisation's [order, theta] arrays.
    K.k rule: past the fringe-parallel angle the conjugate integer order is "order 1"."""
    rows = []
    for j, t in enumerate(sign * th_l):
        flip = math.cos(t * deg + case["psi"]) < 0
        use = kc if (flip and kc is not None) else k
        pr, pt = (arr[use, j], other[use, j]) if reflection else (other[use, j], arr[use, j])
        plus = (arr[k, j], other[k, j]) if reflection else (other[k, j], arr[k, j])
        minus = ((arr[kc, j], other[kc, j]) if reflection else (other[kc, j], arr[kc, j])) if kc is not None else ("", "")
        rows.append(("angle", -float(t), case["lam_um"], float(pr), float(pt), int(orders[use]), harmonics, layers,
                     plus[0], plus[1], minus[0], minus[1]))
    return rows


def _preview(rows, case, reflection, pol):
    """Gate L1 preview (the authority is compare.py / summary<SUF>.csv)."""
    ang = np.array([(-r[1], r[3] if reflection else r[4]) for r in rows if r[0] == "angle"])
    ref = referee(case, ang[:, 0], case["lam_nm"], reflection, pol=pol)
    refl = referee(case, ang[:, 0], case["lam_nm"], reflection, "linear", pol=pol)
    m = np.isfinite(ref)
    print(f"   {pol} angle sweep vs referee: exact rms {np.sqrt(np.mean((ang[m, 1] - ref[m]) ** 2)):.2e}, "
          f"eps-linear rms {np.sqrt(np.mean((ang[m, 1] - refl[m]) ** 2)):.2e}   (gate L1: <= 1e-3 against one of them)", flush=True)


def sweep(fdtd, case, name, harmonics, layers, full):
    """Both sweeps; returns (rows_TE, rows_TM) from the SAME solves (2026-10-16: the p attributes
    are read beside the s ones; the order selection and the theta-sign fit are made on TE, where
    the referee has been confirmed, and applied unchanged to TM)."""
    reflection = MODEL[name][2] == "R"
    # ---- angle sweep at the design wavelength
    th = case["theta_deg"]
    configure(fdtd, case, name, th, [case["lam_um"] * 1e-6], harmonics, layers, full)
    res = solve(fdtd, case, full)
    arr = orders_array(res, reflection)[:, 0, :]
    other = orders_array(res, not reflection)[:, 0, :]
    arr_p = orders_array(res, reflection, "TM")[:, 0, :]
    other_p = orders_array(res, not reflection, "TM")[:, 0, :]
    th_l = axis(res, P["res_theta"])
    orders = np.asarray(res[P["res_orders"]]).ravel()
    sign, k = choose_sign(th_l, arr, orders, case, case["lam_nm"], reflection)
    kc = conj_index(orders, k)
    rows = _rows_from(arr, other, reflection, orders, k, kc, sign, th_l, case, harmonics, layers)
    rows_p = _rows_from(arr_p, other_p, reflection, orders, k, kc, sign, th_l, case, harmonics, layers)
    # ---- wavelength sweep at the Bragg angle
    lams = case["lam_sweep_nm"] * 1e-9
    configure(fdtd, case, name, [case["theta_bragg_deg"] * sign], lams, harmonics, layers, full)
    res = solve(fdtd, case, full)
    arr = orders_array(res, reflection)[:, :, 0]
    other = orders_array(res, not reflection)[:, :, 0]
    arr_p = orders_array(res, reflection, "TM")[:, :, 0]
    other_p = orders_array(res, not reflection, "TM")[:, :, 0]
    lam_l = axis(res, P["res_lambda"], fallback="f")
    # grating_power only carries the orders that PROPAGATE somewhere in the sweep, so this solve's
    # order axis need not be the angle sweep's.  Find the Bragg integer order in it by value.
    orders_w = np.asarray(res[P["res_orders"]]).ravel()
    hit = np.flatnonzero(orders_w == orders[k])
    if hit.size != 1:
        raise RuntimeError(f"order {orders[k]:.0f} is not in the wavelength sweep's order axis {orders_w.tolist()}")
    kw = int(hit[0])
    for j, l in enumerate(lam_l):
        pr, pt = (arr[kw, j], other[kw, j]) if reflection else (other[kw, j], arr[kw, j])
        rows.append(("wavelength", -float(case["theta_bragg_deg"]), float(l) * 1e6, float(pr), float(pt), int(orders[k]),
                     harmonics, layers, "", "", "", ""))
        pr, pt = (arr_p[kw, j], other_p[kw, j]) if reflection else (other_p[kw, j], arr_p[kw, j])
        rows_p.append(("wavelength", -float(case["theta_bragg_deg"]), float(l) * 1e6, float(pr), float(pt), int(orders[k]),
                       harmonics, layers, "", "", "", ""))
    _preview(rows, case, reflection, "TE")
    _preview(rows_p, case, reflection, "TM")
    return rows, rows_p


def flag(argv, key, default, cast):
    return cast(argv[argv.index(key) + 1]) if key in argv else default


def main(argv):
    name = [a for a in argv if not a.startswith("--") and not a.replace(".", "").replace("-", "").isdigit()][0]
    case = CASES[name]
    harmonics = flag(argv, "--harmonics", None, int)
    layers = flag(argv, "--layers", None, int)
    full = "--full-structure" in argv
    global CHI0
    CHI0 = flag(argv, "--chi0", 0.0, float)
    if CHI0 > 0.0 and "--bragg" not in argv:
        raise SystemExit("--chi0 is an UNTESTED single-point absorption probe: use it with --bragg only")
    lumapi = load_lumapi()
    fdtd = open_model(lumapi, FSP_UNSLANTED if (name == "unslanted_T" and os.path.isfile(FSP_UNSLANTED)) else FSP)
    print("Lumerical version:", fdtd.version(), flush=True)
    rows, rows_tm, failed = [], [], False
    try:
        if "--probe" in argv:
            probe(fdtd, name)
            return
        if "--bragg" in argv:
            bragg_check(fdtd, case, name, harmonics, layers, full)
            return
        if "--converge" in argv:
            converge(fdtd, case, name, full)
            return
        bragg_check(fdtd, case, name, harmonics, layers, full)
        rows, rows_tm = sweep(fdtd, case, name, harmonics, layers, full)
        os.makedirs(os.path.join(ROOT, "lumerical"), exist_ok=True)
        fdtd.save(os.path.join(ROOT, "lumerical", f"rcwa_{name}.fsp"))
    except Exception:
        import traceback
        traceback.print_exc()
        failed = True
    finally:
        if not any(f in argv for f in ("--probe", "--bragg", "--converge")):
            os.makedirs(RES, exist_ok=True)
            head = ["sweep", "source_tilt_x_deg", "lambda_um", "P_entry_side_W", "P_far_side_W", "order",
                    "harmonics", "layers_per_period", "P_entry_plus_W", "P_far_plus_W", "P_entry_minus_W", "P_far_minus_W"]
            for suf, rr in (("", rows), ("_TM", rows_tm)):
                bad = failed or not rr
                out = os.path.join(RES, f"lumerical_rcwa_{name}{suf}{'.FAILED' if bad else ''}.csv")
                with open(out, "w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(head)
                    w.writerows(rr)
                print(f"{len(rr)} rows -> {out}", flush=True)
        try:
            fdtd.close()
        except Exception as e:
            print("close:", e)
        sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1:])
