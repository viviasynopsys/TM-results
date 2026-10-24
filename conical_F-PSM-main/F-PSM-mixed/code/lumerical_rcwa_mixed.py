"""lumerical_rcwa_mixed -- drive Ansys Lumerical's RCWA solver on MIXED phase / absorption volume
gratings at conical incidence (the configs of cases_mixed.py) and write per-config efficiency and
absorption tables that compare_mixed.py scores against the complex-permittivity conical referee
(results/reference_mixed.csv, computed on the benchmark PC).

STATUS 2026-10-23: built on the benchmark PC, which has no Lumerical.  Everything lossless here is
the VERIFIED machinery of the conical leg (lumerical_rcwa.py of 2026-09-18, lumerical_rcwa_conical.py
of 2026-10-23: model, thickness through layer repetition, (theta, phi) table mode, sign fit, powers,
amplitudes).  What is NEW is the complex index.  Two routes:

STATUS 2026-10-23 (Lumerical 2026 R1 build 8.35.4609, solver PC): the `import` route WORKS and is the
only route that does (VERIFIED 2026-10-23, P01 against the exact Fresnel oracle to 1.4e-13, P02 against
the referee to 6.7e-5, then 58 configs of sets A-C against the referee to <= 2.1e-4 rms).  The
`equation` route is REJECTED (VERIFIED 2026-10-23): for every token (1i, i, 1j) and both signs the
rectangle's index equation reads back as the real equation, the imaginary part is not stored.
One change to this file that day: remove_import() calls switchtolayout() before deleting the import
object (the session is in analysis mode after a solve).  Route details below keep their original
wording; the per-route status lines carry the result.

  equation   the rectangle's spatial index equation (the setup script writes "n0+Dn*cos(ARG)") is
             replaced by  n0 + n1 cos(ARG) + (kappa0 + kappa1 cos(ARG + delta)) * 1i.  Lumerical's time
             convention is exp(-i w t) (VERIFIED 2026-10-23 from the amplitudes), so + i kappa absorbs;
             the probe tries + first and flips only if the power sum exceeds 1.  Whether the equation
             parser accepts a complex value is UNKNOWN (a 2025 forum reply says it is real-only); the
             imaginary token is tried as 1i, then i, then 1j.  A parser that silently drops the
             imaginary part returns power sums of exactly 1: the probe detects that and the route fails.
             RESULT 2026-10-23: REJECTED (VERIFIED), the stored equation never carries the imaginary part.
  import     an "import" primitive carrying spatial (n, k) data, importnk2(n, x, y, z) (Lumerical
             script reference: n complex, imaginary part = k, x/y/z uniformly spaced, metres), built
             on a 256 x 2 x 1024 grid over one period of the rectangle, given mesh order 1 so it
             overrides the rectangle, whose own index is set to the plain n0 so that its faces still
             generate the layer interfaces.  RESULT 2026-10-23: the RCWA solver DOES sample the import
             object (VERIFIED 2026-10-23); k sign +1 (n + i kappa absorbs); the mesh-order properties were set without error.

The loss probe (--probe-loss) decides the route, the sign and the token from the data (a uniform
absorbing slab against the exact Fresnel oracle, then the published absorbing point of 2026-10-16),
and records them in results/loss_route.json, which every later run reads.  Nothing is tuned.

Usage (from code/; Windows Python with lumapi reachable, see check_env.py):
    python lumerical_rcwa_mixed.py --probe-loss [--index-route equation|import] [--harmonics 3] [--layers 160]
    python lumerical_rcwa_mixed.py --set A|B|C|D|P [--redo]          # every config of a set; existing good tables are skipped
    python lumerical_rcwa_mixed.py --config <case_id> [--redo]
    python lumerical_rcwa_mixed.py --converge <case_id>              # layers 80/160/320 and harmonics 3/5 at theta0
    python lumerical_rcwa_mixed.py --list
Each config writes results/lumerical_mixed_<case_id>.csv (a failed one *.FAILED.csv; a good table is
never overwritten without --redo).  The model is saved as lumerical/rcwa_mixed_<set>.fsp for the first
config of each set.
"""
import csv
import json
import math
import os
import sys
import traceback

import numpy as np

import _paths  # noqa: F401
import lumerical_rcwa as sept
from lumerical_rcwa import P, MODEL
import lumerical_rcwa_conical as drv
import cases_mixed as cm
from mixed_material import solve_mixed
from fresnel_slab import slab
from rcwa_conical import eta_of

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
LUM = os.path.join(ROOT, "lumerical")
deg = math.pi / 180.0

ROUTES = ("equation", "import")
IMAG_TOKENS = ("1i", "i", "1j")      # equation route only; all three rejected by the solver PC's Lumerical 2026 R1 (VERIFIED 2026-10-23)
IMPORT_NAME = "nkimport"
NX, NZ = 256, 1024                   # import-route grid over one period (x) and one z-period (or the plate, unslanted)
CURRENT = None                       # the config being solved (read by the mock in the dry run)
APPLIED_DELTA_SIGN = +1              # the orientation of delta last written into the model (read by the mock)
AMP_KEYS = drv.AMP_KEYS

HEAD = ["set", "case_id", "family", "grating", "phi_deg", "theta_deg", "lambda_um", "order", "n1", "kap0", "kap1", "delta_deg",
        "d_um", "index_route", "k_sign", "imag_token", "delta_sign", "theta_lum_deg", "phi_lum_deg", "eta_s_lum", "eta_p_lum",
        "power_s_lum", "power_p_lum", "absorb_s_lum", "absorb_p_lum", "harmonics", "layers_per_period", "theta_sign"]
HEAD += [f"{k}_{ri}" for k in AMP_KEYS for ri in ("re", "im")]
COL = {k: i for i, k in enumerate(HEAD)}


# ----------------------------------------------------------------------------- model inputs
def inputs(cfg):
    """(name, case_dict for the solve, case_list for the sign fit).  BENCH gratings use the September
    model inputs; DAVID gratings get recording beams from cases_conical (K-vector residual checked)."""
    c0 = dict(set=cfg["model_set"], case=cfg["grating"], tau=cfg["tau"], k0=2 * math.pi / cfg["lam_um"], n0=cfg["n0"],
              n1=cfg["n1_fit"], K=cfg["K"], Psi=cfg["Psi"], d=cfg["d"], lam_um=cfg["lam_um"], phi_deg=cfg["phi_deg"],
              theta_bragg_deg=cfg["theta0_deg"], thetas_deg=cfg["thetas_deg"])
    case_dict, name = drv.model_inputs([c0])
    case_dict = dict(case_dict)                      # never mutate sept.CASES
    case_dict["d_um"] = cfg["d"]
    case_dict["lam_um"] = cfg["lam_um"]
    fit_dict = dict(case_dict, n1=cfg["n1_fit"])     # the lossless sign fit needs a diffracting grating
    solve_dict = dict(case_dict, n1=cfg["n1"])
    return name, solve_dict, fit_dict, [c0]


# ----------------------------------------------------------------------------- complex index
def cos_argument(expr):
    """The argument of the cos( ... ) in the rectangle's index equation (balanced parentheses)."""
    i = expr.find("cos(")
    if i < 0:
        raise RuntimeError(f"no cos( in the rectangle index equation: {expr!r}")
    depth = 0
    for k in range(i + 3, len(expr)):
        if expr[k] == "(":
            depth += 1
        elif expr[k] == ")":
            depth -= 1
            if depth == 0:
                return expr[i + 4:k]
    raise RuntimeError(f"unbalanced parentheses in the index equation: {expr!r}")


def arg_function(name, case_dict, kz_zero):
    """ARG(x_um, z_um) as the setup script builds it (x, y, z in MICRONS), at full precision.  The
    script rounds the K components with num2str (6 significant figures); check_argument reports the
    difference and refuses above 1e-4 rad."""
    a1, a2, typ = MODEL[name]
    n0, lam_r = case_dict["n0"], case_dict["lam_r_um"]
    if kz_zero:
        dx_um = lam_r / (2.0 * n0 * math.sin(a1 * deg))
        return lambda x, z: 2.0 * math.pi * x / dx_um
    k1 = np.array([math.sin(a1 * deg), math.cos(a1 * deg)])
    k2 = np.array([math.sin(a2 * deg), (-1.0 if typ == "R" else 1.0) * math.cos(a2 * deg)])
    Kv = k1 - k2
    pref = 2.0 * math.pi * n0 / lam_r
    return lambda x, z: pref * (x * Kv[0] + z * Kv[1])


def check_argument(arg_str, fn, x_um, z_um):
    """Evaluate the script's own ARG string (Lumerical syntax is Python syntax here) against fn."""
    env = {"pi": np.pi, "cos": np.cos, "sin": np.sin, "x": x_um, "y": np.zeros_like(x_um), "z": z_um}
    ref = eval(arg_str, {"__builtins__": {}}, env)
    dev = float(np.max(np.abs(ref - fn(x_um, z_um))))
    if dev > 1e-4:
        raise RuntimeError(f"the setup script's cos argument differs from the model inputs by {dev:.2e} rad: {arg_str!r}")
    return dev


def complex_expression(n0, n1, kap0, kap1, delta_rad, arg, ksign, token):
    sign = "+" if ksign > 0 else "-"
    return f"{n0!r}+{n1!r}*cos({arg}){sign}({kap0!r}+{kap1!r}*cos({arg}+{delta_rad!r}))*{token}"


def remove_import(fdtd):
    try:
        n = int(fdtd.getnamednumber(IMPORT_NAME))
    except Exception:
        return
    if n > 0:
        fdtd.switchtolayout()                                  # an earlier solve leaves the session in analysis mode (2026-10-23, Lumerical 2026 R1)
        fdtd.select(IMPORT_NAME)
        fdtd.delete()


def apply_equation_route(fdtd, cfg, ksign, token, delta_sign=+1):
    expr0 = str(fdtd.getnamed(P["rect"], P["rect_index"]))
    arg = cos_argument(expr0)
    expr = complex_expression(cfg["n0"], cfg["n1"], cfg["kap0"], cfg["kap1"], delta_sign * cfg["delta_deg"] * deg, arg, ksign, token)
    fdtd.setnamed(P["model"], P["model_script"], "")          # the solve re-runs the script otherwise (VERIFIED 2026-09-18)
    fdtd.setnamed(P["rect"], P["rect_index"], expr)
    back = str(fdtd.getnamed(P["rect"], P["rect_index"]))
    if back.replace(" ", "") != expr.replace(" ", ""):
        raise RuntimeError(f"the index equation did not take: set {expr!r}, read back {back!r}")
    return dict(route="equation", expression=expr, lossless_expression=expr0)


def apply_import_route(fdtd, cfg, name, case_dict, ksign, kz_zero, delta_sign=+1):
    expr0 = str(fdtd.getnamed(P["rect"], P["rect_index"]))
    arg_str = cos_argument(expr0)
    fn = arg_function(name, case_dict, kz_zero)
    xc, yc, zc = (float(fdtd.getnamed(P["rect"], k)) for k in ("x", "y", "z"))
    dz = float(fdtd.getnamed(P["rect"], P["rect_zspan"]))
    dx = float(fdtd.getnamed(P["rcwa"], P["rcwa_xspan"]))
    x = np.linspace(-dx / 2, dx / 2, NX)
    z = np.linspace(-dz / 2, dz / 2, NZ)
    y = np.array([-1e-6, 1e-6])
    X, Z = np.meshgrid(x, z, indexing="ij")
    Xu, Zu = X * 1e6, Z * 1e6                                  # the script's equation takes microns about the rectangle centre
    dev = check_argument(arg_str, fn, Xu, Zu)
    ARG = fn(Xu, Zu)
    Nidx = cfg["n0"] + cfg["n1"] * np.cos(ARG) + ksign * 1j * (cfg["kap0"] + cfg["kap1"] * np.cos(ARG + delta_sign * cfg["delta_deg"] * deg))
    N3 = np.repeat(Nidx[:, None, :], 2, axis=1)                # (NX, 2, NZ)
    remove_import(fdtd)
    fdtd.addimport()
    fdtd.set("name", IMPORT_NAME)
    fdtd.importnk2(N3, x + xc, y + yc, z + zc)                 # VERIFIED 2026-10-23: signature importnk2(n, x, y, z), metres, accepted
    mesh_note = "mesh order 1 set"
    try:
        fdtd.set("override mesh order from material database", 1)
        fdtd.set("mesh order", 1)
    except Exception as e:
        mesh_note = f"mesh order NOT set: {str(e).splitlines()[0]}"
    got = {k: float(fdtd.getnamed(IMPORT_NAME, k)) for k in ("x", "y", "z", "x span", "z span")}
    if abs(got["x span"] - dx) > 1e-12 or abs(got["z span"] - dz) > 1e-12 or abs(got["x"] - xc) > 1e-12 or abs(got["z"] - zc) > 1e-12:
        raise RuntimeError(f"the import object does not cover the rectangle: {got} against x span {dx}, z span {dz}, centre ({xc}, {zc})")
    fdtd.setnamed(P["model"], P["model_script"], "")
    fdtd.setnamed(P["rect"], P["rect_index"], repr(float(cfg["n0"])))
    return dict(route="import", grid=(NX, 2, NZ), arg_deviation=dev, mesh=mesh_note, lossless_expression=expr0)


def apply_route(fdtd, cfg, name, case_dict, route, ksign, token, delta_sign=+1):
    global APPLIED_DELTA_SIGN
    APPLIED_DELTA_SIGN = delta_sign
    if cfg["kap0"] == 0.0 and cfg["kap1"] == 0.0:
        return dict(route="lossless")
    kz_zero = abs(math.cos(float(case_dict["psi"]))) < 1e-9
    if route == "equation":
        return apply_equation_route(fdtd, cfg, ksign, token, delta_sign)
    if route == "import":
        return apply_import_route(fdtd, cfg, name, case_dict, ksign, kz_zero, delta_sign)
    raise ValueError(f"index route must be one of {ROUTES}, not {route!r}")


# ----------------------------------------------------------------------------- one solve
def solve_rows(fdtd, cfg, name, case_dict, pairs, harmonics, layers, route, ksign, token, delta_sign=+1):
    """Configure the lossless model (verified path), apply the complex index, solve, read everything."""
    global CURRENT
    CURRENT = cfg
    remove_import(fdtd)
    drv.configure_conical(fdtd, case_dict, name, pairs, cfg["lam_um"] * 1e-6, harmonics, layers)
    info = apply_route(fdtd, cfg, name, case_dict, route, ksign, token, delta_sign)
    res = sept.solve(fdtd, case_dict, False)
    orders = np.asarray(res[P["res_orders"]]).ravel()
    r_s, t_s = drv.read_powers(res, "s")
    r_p, t_p = drv.read_powers(res, "p")
    er_s, et_s = sept.total_energy(fdtd, True, "TE")
    er_p, et_p = sept.total_energy(fdtd, True, "TM")
    pow_s = np.asarray(er_s).ravel() + np.asarray(et_s).ravel()
    pow_p = np.asarray(er_p).ravel() + np.asarray(et_p).ravel()
    amps, amp_orders = drv.read_amplitudes(fdtd)
    return dict(orders=orders, r_s=r_s, t_s=t_s, r_p=r_p, t_p=t_p, pow_s=pow_s, pow_p=pow_p, amps=amps, amp_orders=amp_orders, info=info)


_SIGN_CACHE = {}


def sign_fit(fdtd, cfg, name, fit_dict, case_list, harmonics, layers):
    """Lumerical's theta sign and the signal order, fitted once per grating geometry on the lossless
    planar replay (lumerical_rcwa_conical.planar_sign_fit, VERIFIED 2026-10-23)."""
    global CURRENT
    CURRENT = cfg                                   # the config in hand, also during the fit (the dry run's mock reads it)
    key = (name, round(cfg["K"], 9), round(cfg["Psi"], 9))
    if key not in _SIGN_CACHE:
        remove_import(fdtd)
        _SIGN_CACHE[key] = drv.planar_sign_fit(fdtd, case_list, fit_dict, name, harmonics, layers)
    return _SIGN_CACHE[key]


_DELTA_CACHE = {}
_DELTA_LOG = []


def needs_delta(cfg):
    """The sign of delta matters only when both harmonics are present and delta is not 0 or 180 deg:
    with n1 = 0 the relative phase is a shift of the whole profile (|B| = |C|), and at 0 / 180 deg
    B and C have equal magnitude.  (Mixed report Sec. 9.A: C != B* in general.)"""
    return cfg["kap1"] != 0.0 and cfg["n1"] != 0.0 and abs(math.fmod(cfg["delta_deg"], 180.0)) > 1e-9


def delta_orientation(fdtd, cfg, harmonics, layers, route, ksign, token):
    """Which way round Lumerical reads the relative phase delta for this grating geometry, from the
    data.  The lossless planar sign fit cannot decide it: the mirror x' -> -x' of the slab maps
    (theta sign, order, delta) -> (-sign, -order, -delta), and for an unslanted grating the two
    images fit the lossless curve equally well (a tie the dry run exposed); for a slanted one the
    fit is unique but which exponential, e^{+i ARG} or e^{-i ARG}, feeds Lumerical's order +1 is
    still a convention of the solver.  So, once per geometry, the ONE-SIDED grating of the author's
    Sec. 9.D is solved at theta0 (n1 = kappa1 = 0.012, kappa0 = 0.015, d = 2 um, delta = 270 deg):
    with the orientation the referee uses, order +1 is dark (exactly 0, the A10/A20/A30 cases);
    with the opposite orientation it carries the delta = 90 deg efficiency (0.02 .. 0.04).  The
    ratio is 1e4 or more, so the decision is unambiguous; anything in between is refused."""
    key = (cfg["grating"], round(cfg["K"], 9), round(cfg["Psi"], 9))
    if key in _DELTA_CACHE:
        return _DELTA_CACHE[key]
    pc = dict(cfg, case_id=cfg["case_id"] + "_DELTAPROBE", n1=0.012, n1_fit=0.012, kap1=0.012, kap0=0.015, delta_deg=270.0, d=2.0,
              thetas_deg=np.array([cfg["theta0_deg"]]), bound_E19="", needs_sign_fit=True)
    name, solve_dict, fit_dict, case_list = inputs(pc)
    sign, k, orders_fit = sign_fit(fdtd, pc, name, fit_dict, case_list, harmonics, layers)
    pairs = drv.lum_pairs(pc["thetas_deg"], pc["phi_deg"], sign)
    S = solve_rows(fdtd, pc, name, solve_dict, pairs, harmonics, layers, route, ksign, token, +1)
    kk = int(np.flatnonzero(S["orders"] == orders_fit[k])[0])
    got = float((S["r_s"] if pc["tau"] == -1 else S["t_s"])[kk, 0])
    k0 = 2 * math.pi / pc["lam_um"]
    th = pc["theta0_deg"] * deg
    dark = solve_mixed(k0, pc["n0"], pc["K"], pc["Psi"], th, pc["phi_deg"] * deg, pc["d"], pc["n1"], pc["kap0"], pc["kap1"], 270.0, M=10)
    lit = solve_mixed(k0, pc["n0"], pc["K"], pc["Psi"], th, pc["phi_deg"] * deg, pc["d"], pc["n1"], pc["kap0"], pc["kap1"], 90.0, M=10)
    e_dark = eta_of(dark["J_T"] if pc["tau"] == +1 else dark["J_R"], [1, 0])
    e_lit = eta_of(lit["J_T"] if pc["tau"] == +1 else lit["J_R"], [1, 0])
    if got < 0.1 * e_lit:
        dsign = +1
    elif abs(got - e_lit) < 0.1 * e_lit:
        dsign = -1
    else:
        raise RuntimeError(f"delta orientation probe on {cfg['grating']}: Lumerical gives {got:.6f} for the one-sided grating, "
                           f"referee {e_dark:.2e} with +delta and {e_lit:.6f} with -delta: neither; stop and report")
    _DELTA_CACHE[key] = dsign
    entry = dict(grating=cfg["grating"], phi_deg=cfg["phi_deg"], theta0_deg=cfg["theta0_deg"], probe_case=pc["case_id"],
                 lumerical_eta_s=got, referee_eta_s_plus_delta=e_dark, referee_eta_s_minus_delta=e_lit, delta_sign=dsign,
                 theta_sign=int(sign), signal_order=int(orders_fit[k]), status="VERIFIED from the data")
    _DELTA_LOG.append(entry)
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "delta_orientation.json"), "w") as f:
        json.dump(_DELTA_LOG, f, indent=1)
    print(f"   delta orientation for {cfg['grating']}: one-sided probe gives eta_s {got:.6f} (referee: {e_dark:.1e} with +delta, "
          f"{e_lit:.6f} with -delta) -> delta sign {dsign:+d}", flush=True)
    return dsign


# ----------------------------------------------------------------------------- reference preview
_REF = None


def reference():
    global _REF
    if _REF is None:
        _REF = {}
        path = os.path.join(RES, "reference_mixed.csv")
        if os.path.isfile(path):
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    _REF[(r["case_id"], round(float(r["theta_deg"]), 6))] = r
    return _REF


def preview(cfg, thetas, got_s, got_p, abs_s, abs_p):
    ref = reference()
    rs, rp, ra_s, ra_p = [], [], [], []
    k0 = 2 * math.pi / cfg["lam_um"]
    for th in thetas:
        r = ref.get((cfg["case_id"], round(float(th), 6)))
        if r is not None and cfg["set"] != "P":
            rs.append(float(r["eta_s"])); rp.append(float(r["eta_p"]))
            ra_s.append(float(r["absorb_s"])); ra_p.append(float(r["absorb_p"]))
        else:
            R = solve_mixed(k0, cfg["n0"], cfg["K"], cfg["Psi"], th * deg, cfg["phi_deg"] * deg, cfg["d"], cfg["n1"], cfg["kap0"],
                            cfg["kap1"], cfg["delta_deg"], M=10)
            J = R["J_T"] if cfg["tau"] == +1 else R["J_R"]
            if cfg["set"] == "P" and cfg["needs_sign_fit"] is False:
                rs.append(float(R["DE_T"][10, 0])); rp.append(float(R["DE_T"][10, 1]))
            else:
                rs.append(eta_of(J, [1, 0])); rp.append(eta_of(J, [0, 1]))
            ra_s.append(R["absorb_s"]); ra_p.append(R["absorb_p"])
    rms = lambda a, b: float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))
    print(f"   {cfg['case_id']}: s rms vs referee {rms(got_s, rs):.2e}, p rms {rms(got_p, rp):.2e}; absorbed-fraction rms s {rms(abs_s, ra_s):.2e} "
          f"p {rms(abs_p, ra_p):.2e}; power sums s {min(1 - np.asarray(abs_s)):.6f}..{max(1 - np.asarray(abs_s)):.6f}   "
          f"(gates: rms <= 1e-3, sums < 1)", flush=True)


# ----------------------------------------------------------------------------- the loss probe
def probe_loss(fdtd, route, harmonics, layers):
    log = [f"Lumerical loss probe, route {route!r}, {harmonics} harmonics, {layers} layers per period; Lumerical {fdtd.version()}\n"]

    def say(s):
        log.append(s + "\n")
        print(s, flush=True)

    cfg = cm.config("P01_lumerical_R_n000_k00000_02000_d000_phi150")
    name, solve_dict, fit_dict, case_list = inputs(cfg)
    pairs = drv.lum_pairs(cfg["thetas_deg"], cfg["phi_deg"], +1.0)
    k0 = 2 * math.pi / cfg["lam_um"]
    oracle = [slab(k0, cfg["n0"], cfg["kap0"], th * deg, cfg["phi_deg"] * deg, cfg["d"]) for th in cfg["thetas_deg"]]
    chosen = None
    tokens = IMAG_TOKENS if route == "equation" else ("n+ik array",)
    for token in tokens:
        for ksign in (+1, -1):
            say(f"\nP01 uniform slab n0 - i {cfg['kap0']} (d {cfg['d']} um): route {route}, token {token}, Lumerical k sign {ksign:+d}")
            try:
                S = solve_rows(fdtd, cfg, name, solve_dict, pairs, harmonics, layers, route, ksign, token)
            except Exception as e:
                say(f"   FAILED: {type(e).__name__}: {str(e).splitlines()[0] if str(e) else ''}")
                continue
            k0i = int(np.flatnonzero(S["orders"] == 0)[0])
            sums = np.concatenate([S["pow_s"], S["pow_p"]])
            say(f"   power sums s {S['pow_s']} p {S['pow_p']}")
            for j, th in enumerate(cfg["thetas_deg"]):
                say(f"   theta {th:5.1f}: T0 s {S['t_s'][k0i, j]:.6f} (oracle {oracle[j]['s'][0]:.6f})  p {S['t_p'][k0i, j]:.6f} ({oracle[j]['p'][0]:.6f})   "
                    f"R0 s {S['r_s'][k0i, j]:.6f} ({oracle[j]['s'][1]:.6f})  p {S['r_p'][k0i, j]:.6f} ({oracle[j]['p'][1]:.6f})")
            if np.max(np.abs(sums - 1.0)) < 1e-6:
                say("   the power sums are exactly 1: the imaginary part was IGNORED on this route/token")
                break                                        # the sign cannot help; try the next token
            if np.min(sums) > 1.0 + 1e-6:
                say("   the power sums exceed 1: this sign is GAIN in Lumerical's convention; trying the other sign")
                continue
            dev = max(abs(S["t_s"][k0i, j] - oracle[j]["s"][0]) for j in range(len(cfg["thetas_deg"])))
            dev = max(dev, max(abs(S["t_p"][k0i, j] - oracle[j]["p"][0]) for j in range(len(cfg["thetas_deg"]))))
            say(f"   absorbing; max |T0 - Fresnel oracle| = {dev:.2e} (gate 1e-3)")
            if dev <= 1e-3:
                chosen = dict(route=route, k_sign=ksign, imag_token=token, p01_max_dev=dev)
                break
            say("   the loss is applied but does not match the oracle: route rejected")
        if chosen:
            break
    if chosen is None:
        with open(os.path.join(RES, f"lumerical_mixed_probe_loss_{route}.FAILED.txt"), "w", encoding="utf-8") as f:
            f.writelines(log)
        raise SystemExit(f"loss probe FAILED on route {route!r}: no (sign, token) applies the imaginary index and matches the Fresnel slab. "
                         f"See results/lumerical_mixed_probe_loss_{route}.FAILED.txt; try the other route, then stop and report.")
    say(f"\nP01 established (VERIFIED from the data): {chosen}")

    # P02: the published absorbing point
    cfg2 = cm.config("P02_lumerical_R_n040_k00000_00100_d000_phi000")
    name, solve_dict, fit_dict, case_list = inputs(cfg2)
    sign, k, orders = sign_fit(fdtd, cfg2, name, fit_dict, case_list, harmonics, layers)
    pairs = drv.lum_pairs(cfg2["thetas_deg"], cfg2["phi_deg"], sign)
    S = solve_rows(fdtd, cfg2, name, solve_dict, pairs, harmonics, layers, chosen["route"], chosen["k_sign"], chosen["imag_token"])
    kk = int(np.flatnonzero(S["orders"] == orders[k])[0])
    te, tm = S["r_s"][kk, 0], S["r_p"][kk, 0]
    say(f"\nP02 lumerical_R theta 0, kappa0 1e-3: R+1 TE {te:.6f} TM {tm:.6f}; power sums s {S['pow_s'][0]:.6f} p {S['pow_p'][0]:.6f}"
        f"   (published referee: 0.735964 / 0.668751, sums 0.915638 / 0.909354; gate 2e-3 with the 160-layer staircase)")
    worst = max(abs(te - 0.735964), abs(tm - 0.668751), abs(S["pow_s"][0] - 0.915638), abs(S["pow_p"][0] - 0.909354))
    ok2 = worst <= 2e-3
    say(f"P02 max |deviation| {worst:.2e}: {'PASS' if ok2 else 'FAIL'}")
    chosen.update(p02_max_dev=float(worst), p02_pass=bool(ok2), harmonics=harmonics, layers=layers,
                  lumerical_version=str(fdtd.version()), date="2026-10-23", status="VERIFIED from the data on this PC")
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "lumerical_mixed_probe_loss.txt"), "w", encoding="utf-8") as f:
        f.writelines(log)
    with open(os.path.join(RES, "loss_route.json"), "w") as f:
        json.dump(chosen, f, indent=1)
    print(f"-> results/lumerical_mixed_probe_loss.txt, results/loss_route.json", flush=True)
    if not ok2:
        raise SystemExit("P02 failed: the loss route applies absorption but the published point is not reproduced. Stop and report.")


# ----------------------------------------------------------------------------- sweeps
def run_config(fdtd, cfg, harmonics, layers, route, ksign, token, redo=False, save_model=False):
    global CURRENT
    CURRENT = cfg
    out = os.path.join(RES, f"lumerical_mixed_{cfg['case_id']}.csv")
    if os.path.isfile(out) and not redo:
        print(f"   {cfg['case_id']}: table exists, skipped (use --redo)", flush=True)
        return
    rows, failed = [], False
    try:
        name, solve_dict, fit_dict, case_list = inputs(cfg)
        if cfg["needs_sign_fit"]:
            sign, k, orders_fit = sign_fit(fdtd, cfg, name, fit_dict, case_list, harmonics, layers)
            want_order = int(orders_fit[k])
        else:
            sign, want_order = +1.0, 0
        dsign = delta_orientation(fdtd, cfg, harmonics, layers, route, ksign, token) if needs_delta(cfg) else +1
        pairs = drv.lum_pairs(cfg["thetas_deg"], cfg["phi_deg"], sign)
        S = solve_rows(fdtd, cfg, name, solve_dict, pairs, harmonics, layers, route, ksign, token, dsign)
        hit = np.flatnonzero(S["orders"] == want_order)
        if hit.size != 1:
            raise RuntimeError(f"order {want_order:+d} not in this solve's order axis {S['orders'].tolist()}")
        kk = int(hit[0])
        ka = None
        if S["amps"] is not None:
            ha = np.flatnonzero(S["amp_orders"] == want_order)
            ka = int(ha[0]) if ha.size == 1 else None
        reflect = cfg["tau"] == -1 and cfg["needs_sign_fit"]
        for j, th in enumerate(cfg["thetas_deg"]):
            es = (S["r_s"] if reflect else S["t_s"])[kk, j]
            ep = (S["r_p"] if reflect else S["t_p"])[kk, j]
            row = [cfg["set"], cfg["case_id"], cfg["family"], cfg["grating"], float(cfg["phi_deg"]), float(th), cfg["lam_um"], want_order,
                   cfg["n1"], cfg["kap0"], cfg["kap1"], cfg["delta_deg"], cfg["d"], S["info"]["route"], ksign, token, dsign,
                   float(pairs[j, 0]), float(pairs[j, 1]), float(es), float(ep), float(S["pow_s"][j]), float(S["pow_p"][j]),
                   float(1.0 - S["pow_s"][j]), float(1.0 - S["pow_p"][j]), harmonics, layers, int(sign)]
            if S["amps"] is not None and ka is not None:
                for key in AMP_KEYS:
                    v = S["amps"][key][ka, j]
                    row += [float(np.real(v)), float(np.imag(v))]
            else:
                row += [""] * (2 * len(AMP_KEYS))
            rows.append(row)
        preview(cfg, cfg["thetas_deg"], [r[COL["eta_s_lum"]] for r in rows], [r[COL["eta_p_lum"]] for r in rows],
                [r[COL["absorb_s_lum"]] for r in rows], [r[COL["absorb_p_lum"]] for r in rows])
        if save_model:
            os.makedirs(LUM, exist_ok=True)
            fdtd.save(os.path.join(LUM, f"rcwa_mixed_{cfg['set']}.fsp"))
    except Exception:
        traceback.print_exc()
        failed = True
    finally:
        os.makedirs(RES, exist_ok=True)
        bad = failed or not rows
        path = out if not bad else out.replace(".csv", ".FAILED.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(HEAD)
            w.writerows(rows)
        print(f"{len(rows)} rows -> {path}", flush=True)


def converge(fdtd, cfg, route, ksign, token):
    global CURRENT
    CURRENT = cfg
    name, solve_dict, fit_dict, case_list = inputs(cfg)
    lines = [f"convergence ladder, {cfg['case_id']}, theta0 {cfg['theta0_deg']} deg, phi {cfg['phi_deg']} deg\n"]
    k0 = 2 * math.pi / cfg["lam_um"]
    R = solve_mixed(k0, cfg["n0"], cfg["K"], cfg["Psi"], cfg["theta0_deg"] * deg, cfg["phi_deg"] * deg, cfg["d"], cfg["n1"], cfg["kap0"],
                    cfg["kap1"], cfg["delta_deg"], M=10)
    J = R["J_T"] if cfg["tau"] == +1 else R["J_R"]
    lines.append(f"referee (21 orders): eta_s {eta_of(J, [1, 0]):.7f} eta_p {eta_of(J, [0, 1]):.7f} absorbed s {R['absorb_s']:.7f} p {R['absorb_p']:.7f}\n")
    for h, l in ((3, 80), (3, 160), (3, 320), (5, 160)):
        try:
            sign, k, orders_fit = sign_fit(fdtd, cfg, name, fit_dict, case_list, h, l)
            dsign = delta_orientation(fdtd, cfg, h, l, route, ksign, token) if needs_delta(cfg) else +1
            pairs = drv.lum_pairs([cfg["theta0_deg"]], cfg["phi_deg"], sign)
            S = solve_rows(fdtd, cfg, name, solve_dict, pairs, h, l, route, ksign, token, dsign)
            kk = int(np.flatnonzero(S["orders"] == orders_fit[k])[0])
            es = (S["r_s"] if cfg["tau"] == -1 else S["t_s"])[kk, 0]
            ep = (S["r_p"] if cfg["tau"] == -1 else S["t_p"])[kk, 0]
            lines.append(f"H{h} L{l}: eta_s {es:.7f} eta_p {ep:.7f} absorbed s {1 - S['pow_s'][0]:.7f} p {1 - S['pow_p'][0]:.7f}\n")
        except Exception as e:
            lines.append(f"H{h} L{l}: FAILED {str(e).splitlines()[0]}\n")
        print(lines[-1], end="", flush=True)
    path = os.path.join(RES, f"lumerical_mixed_converge_{cfg['case_id']}.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)
    print(f"-> {path}", flush=True)


def route_settings(argv):
    route = sept.flag(argv, "--index-route", None, str)
    ksign = sept.flag(argv, "--k-sign", None, int)
    token = sept.flag(argv, "--imag-token", None, str)
    path = os.path.join(RES, "loss_route.json")
    if os.path.isfile(path):
        with open(path) as f:
            saved = json.load(f)
        route = route or saved["route"]
        ksign = saved["k_sign"] if ksign is None else ksign
        token = token or saved["imag_token"]
    if route is None or ksign is None or token is None:
        raise SystemExit("no loss route established: run  python lumerical_rcwa_mixed.py --probe-loss  first (it writes results/loss_route.json)")
    return route, ksign, token


def main(argv):
    harmonics = sept.flag(argv, "--harmonics", 3, int)
    layers = sept.flag(argv, "--layers", 160, int)
    if "--list" in argv:
        cm.print_table(cm.all_configs())
        return
    lumapi = sept.load_lumapi()
    fdtd = sept.open_model(lumapi, sept.FSP)
    print("Lumerical version:", fdtd.version(), flush=True)
    try:
        if "--probe-loss" in argv:
            probe_loss(fdtd, sept.flag(argv, "--index-route", "equation", str), harmonics, layers)
            return
        route, ksign, token = route_settings(argv)
        print(f"complex index route {route!r}, Lumerical k sign {ksign:+d}, token {token!r} (from results/loss_route.json or flags)", flush=True)
        if "--converge" in argv:
            converge(fdtd, cm.config(sept.flag(argv, "--converge", None, str)), route, ksign, token)
            return
        if "--config" in argv:
            cfgs = [cm.config(sept.flag(argv, "--config", None, str))]
        elif "--set" in argv:
            cfgs = cm.by_set(sept.flag(argv, "--set", None, str).upper())
        else:
            raise SystemExit("give --probe-loss, --set <letter>, --config <case_id>, --converge <case_id> or --list")
        for i, cfg in enumerate(cfgs):
            print(f"\n=== {cfg['case_id']} ({i + 1} of {len(cfgs)}) ===", flush=True)
            run_config(fdtd, cfg, harmonics, layers, route, ksign, token, redo="--redo" in argv, save_model=(i == 0))
    finally:
        try:
            fdtd.close()
        except Exception as e:
            print("close:", e)


if __name__ == "__main__":
    main(sys.argv[1:])
