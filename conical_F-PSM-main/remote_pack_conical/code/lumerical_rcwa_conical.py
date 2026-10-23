"""Drive the Ansys Lumerical RCWA solver at CONICAL incidence for the cases of cases_conical.py
and write per-azimuth efficiency tables that compare_conical.py scores against the conical referee
(results/reference_conical.csv, computed on the benchmark PC with rcwa_conical.py).

Built 2026-10-23 on top of the September/October driver ``lumerical_rcwa.py`` (VERIFIED at the
solver 2026-09-18; its NAME MAP ``P`` and model handling are imported, not copied).  Conical
incidence through the "table" mode, the (theta, phi) map, the powers and the s/p conversion are
VERIFIED 2026-10-23 at the solver (Lumerical 8.35.4609, v261, 2-D Y-normal region, no 3-D fallback
needed): see results/lumerical_conical_verdict.md.  The names it relies on were all seen
in the 2026-09-18 probe dump (reference/lumerical_probe_2026-09-18.txt):

  * the RCWA object's "incident angle" = "table" mode with an N x 2 "incident angle table" whose
    columns are theta and phi in degrees, and the read-only "incident angles actual";
  * result "grating_power" with Rs_power, Ts_power, Rp_power, Tp_power, shaped
    [lambda, (theta,phi) row, n, m, forward|backward];
  * result "grating_characterization" with the complex amplitudes Rss, Rsp, Rps, Rpp, Tss, Tsp,
    Tps, Tpp, same shape.  VERIFIED 2026-10-23: the amplitudes are power-normalised (|amplitude|^2 is
    the order's power) and the order is "in-first" (Rsp = incident s, outgoing p); the s-p relative
    phase of the two sources is NOT reproduced by a +-1 sign (see the verdict) and stays UNTESTED.
    The driver records them raw and compare_conical.py tests, on the planar rows,
    whether |Tss|^2 (or Rss) equals Ts_power (Rs_power) directly or after a flux factor.

Geometry map.  cases_conical.py states every case in the report's convention (theta >= 0 from the
normal, phi from +x' toward +z, slab normal y).  Lumerical's model has x tangential, z normal,
phi from +x toward +y.  The sign of Lumerical's theta relative to the benchmark's was fitted in
September (+1 on every case: Lumerical's theta IS cases.py's theta); it is fitted again here on the
planar rows of each run (phi = 0 and 180 deg are the planar replays at +theta and -theta), and the
same map is applied to phi:  sign +1: (theta_lum, phi_lum) = (theta, phi);  sign -1: (theta, 180 - phi).
Efficiencies are invariant under phi -> -phi (mirror in the plane of incidence), so no further
freedom exists.  The signal order is chosen as in September: the x order (+1 or -1) whose planar
TE curve fits the referee.

Usage (from code/):
    python lumerical_rcwa_conical.py --probe
    python lumerical_rcwa_conical.py --bragg <case>       # planar sign fit + one solve at every azimuth's Bragg point
    python lumerical_rcwa_conical.py <case> [--harmonics 3] [--layers 160]
                                                          # every azimuth sweep -> results/lumerical_conical_<case>.csv
    <case>: lumerical_R, lumerical_T, zemax_kba_R, unslanted_T (BENCH), Slanted_T, Slanted_R, Unslanted_T (DAVID)
A failed run writes *.FAILED.csv and never overwrites a good table.
"""
import csv
import math
import os
import sys

import numpy as np

import lumerical_rcwa as sept                       # the verified September driver
from lumerical_rcwa import P, MODEL
import cases_conical as cc
from rcwa_conical import solve as referee_solve, eta_of

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
deg = np.pi / 180.0

AMP_KEYS = ("Rss", "Rsp", "Rps", "Rpp", "Tss", "Tsp", "Tps", "Tpp")   # names VERIFIED 2026-09-18 and 2026-10-23; power normalisation VERIFIED 2026-10-23
RESULT_CHAR = "grating_characterization"                              # VERIFIED 2026-09-18


def all_cases():
    out = {}
    for c in cc.bench_cases() + cc.david_cases():
        out.setdefault(c["case"], []).append(c)
    return out


def model_inputs(case_list):
    """What set_model needs: a cases.py-like dict and a MODEL entry.  BENCH cases already have
    both; DAVID cases get recording beams from cases_conical.recording_beams."""
    c0 = case_list[0]
    name = c0["case"]
    if c0["set"] == "BENCH":
        return sept.CASES[name], name
    beta = c0["k0"] * c0["n0"]
    a1, a2, typ = cc.recording_beams(c0["K"], c0["Psi"], beta)
    res = cc.check_recording_beams(c0["K"], c0["Psi"], beta, a1, a2, typ)
    if res > 1e-9:
        raise RuntimeError(f"recording beams do not reproduce the grating vector (residual {res:.1e})")
    MODEL[name] = (a1, a2, typ)
    fake = dict(n0=c0["n0"], n1=c0["n1"], lam_r_um=c0["lam_um"], lam_um=c0["lam_um"], d_um=c0["d"],
                K=c0["K"] * 1e-3, psi=c0["Psi"], theta_bragg_deg=c0["theta_bragg_deg"])
    print(f"   DAVID grating {name}: recording beams {a1:.4f} / {a2:.4f} deg, type {typ}, w_rec {c0['lam_um']} um, "
          f"K-vector residual {res:.1e}", flush=True)
    return fake, name


def set_excitation_conical(fdtd, pairs_deg, lam_m):
    """pairs_deg: (N, 2) array of (theta_lum, phi_lum) in degrees; one wavelength."""
    pairs = np.asarray(pairs_deg, float)
    fdtd.setnamed(P["rcwa"], P["theta_mode"], "table")
    fdtd.setnamed(P["rcwa"], P["theta_table"], pairs)
    got = np.asarray(fdtd.getnamed(P["rcwa"], P["theta_actual"]))
    if got.shape != pairs.shape or not np.allclose(got, pairs, atol=1e-9):
        raise RuntimeError(f"asked {pairs.shape[0]} (theta, phi) rows, solver will use {got.shape}: {got[:4]} ...")
    fdtd.setnamed(P["rcwa"], P["use_wl_spacing"], 1)
    fdtd.setnamed(P["rcwa"], P["wl_center"], float(lam_m))
    fdtd.setnamed(P["rcwa"], P["wl_span"], 0.0)
    fdtd.setnamed(P["rcwa"], P["freq_points"], 1)


def configure_conical(fdtd, case_dict, name, pairs_deg, lam_m, harmonics, layers):
    sept.require_names()
    # K_z = 0 (Psi = 90 deg): the setup script's dz is infinite, so use the September bypass.  It reads
    # the recording angle under the key "unslanted_T"; a DAVID family supplies its own beams under its
    # own name.  Added 2026-10-23 on the work PC: the name test alone sent DAVID "Unslanted_T" down the
    # setup-script path, which fails with "asked 160 layers, solver made 41 interfaces".
    kz_zero = abs(math.cos(float(case_dict["psi"]))) < 1e-9
    if name == "unslanted_T" or kz_zero:
        if name != "unslanted_T":
            sept.MODEL["unslanted_T"] = sept.MODEL[name]
        sept.configure_unslanted(fdtd, case_dict, [0.0], [lam_m], harmonics, layers)
        set_excitation_conical(fdtd, pairs_deg, lam_m)
        return
    sept.set_model(fdtd, case_dict, name)
    dz = float(fdtd.getnamed(P["rect"], P["rect_zspan"]))
    sept.set_thickness(fdtd, case_dict, dz, False)
    set_excitation_conical(fdtd, pairs_deg, lam_m)
    sept.set_numerics(fdtd, harmonics, layers, 1.0)
    d_sim = sept.read_thickness(fdtd, False)
    print(f"   z-period {dz * 1e6:.6f} um; simulated thickness {d_sim * 1e6:.4f} um (asked {case_dict['d_um']:.4f})", flush=True)
    if abs(d_sim - case_dict["d_um"] * 1e-6) > 5e-9:
        raise RuntimeError(f"thickness mismatch: solver will simulate {d_sim * 1e6:.4f} um, asked {case_dict['d_um']} um")


def read_powers(res, pol):
    """[order, row] efficiencies for forward incidence, entry side and far side, for incident
    polarisation s (TE) or p (TM)."""
    kR, kT = ("res_R", "res_T") if pol == "s" else ("res_Rp", "res_Tp")
    fwd = int(np.argmax(np.asarray(res[P["res_fwd"]]).ravel()))
    r = np.real(np.asarray(res[P[kR]]))[0, :, :, 0, fwd].T
    t = np.real(np.asarray(res[P[kT]]))[0, :, :, 0, fwd].T
    return r, t


def read_amplitudes(fdtd):
    """Raw complex order amplitudes [key][order, row], forward incidence; None if unavailable."""
    try:
        ch = fdtd.getresult(P["rcwa"], RESULT_CHAR)
    except Exception as e:
        print(f"   {RESULT_CHAR} not readable: {str(e).splitlines()[0]}", flush=True)
        return None, None
    fwd = int(np.argmax(np.asarray(ch[P["res_fwd"]]).ravel()))
    out = {}
    for k in AMP_KEYS:
        if k not in ch:
            print(f"   {RESULT_CHAR} has no {k!r}; amplitudes skipped", flush=True)
            return None, None
        out[k] = np.asarray(ch[k])[0, :, :, 0, fwd].T
    orders = np.asarray(ch[P["res_orders"]]).ravel()
    return out, orders


def lum_pairs(thetas_deg, phi_deg, sign):
    th = np.asarray(thetas_deg, float)
    ph = np.full_like(th, phi_deg if sign > 0 else 180.0 - phi_deg)
    return np.column_stack([th, ph])


def planar_sign_fit(fdtd, case_list, case_dict, name, harmonics, layers):
    """Solve the planar replays (phi = 0 -> +theta, phi = 180 -> -theta in Lumerical's own axes) at
    a few angles about the case's planar Bragg point and fit the theta sign and the signal order
    against the conical referee at phi = 0 / 180.  Returns (sign, order index, orders)."""
    c0 = case_list[0]
    k0, n0, n1, K, Psi, d, tau = (c0[k] for k in ("k0", "n0", "n1", "K", "Psi", "d", "tau"))
    if c0["set"] == "BENCH":
        thb = case_dict["theta_bragg_deg"]
    else:
        # DAVID: c0's theta_bragg_deg is the conical Bragg angle; the planar replay must sit on the
        # planar one, signed in Lumerical's axes (phi_lum = 0 -> +theta, phi = 180 -> -theta).
        # Changed 2026-10-23 on the work PC: sampling about the conical angle put the fit off-Bragg.
        for phi_pl, sgn in ((0.0, 1.0), (math.pi, -1.0)):
            t_pl = cc.bragg_theta_with(k0 * n0, K, Psi, phi_pl)
            if t_pl is not None:
                thb = sgn * t_pl / deg
                break
        else:
            raise RuntimeError("no planar Bragg angle exists for the sign fit")
    thetas = np.array([thb - 3.0, thb - 1.5, thb, thb + 1.5, thb + 3.0])
    pairs = np.column_stack([thetas, np.zeros_like(thetas)])      # Lumerical phi = 0, signed theta
    configure_conical(fdtd, case_dict, name, pairs, c0["lam_um"] * 1e-6, harmonics, layers)
    res = sept.solve(fdtd, case_dict, False)
    orders = np.asarray(res[P["res_orders"]]).ravel()
    r_s, t_s = read_powers(res, "s")
    best = None
    for sign in (+1.0, -1.0):
        ref = []
        for th in sign * thetas:
            theta, phi = (abs(th) * deg, 0.0 if th >= 0 else np.pi)
            R = referee_solve(k0, n0, n1, K, Psi, theta, phi, d, M=10)
            J = R["J_T"] if tau == +1 else R["J_R"]
            ref.append(eta_of(J, [1, 0]))
        ref = np.array(ref)
        for k in range(orders.size):
            if abs(orders[k]) != 1:
                continue
            e = (r_s if tau == -1 else t_s)[k]
            rms = float(np.sqrt(np.mean((e - ref) ** 2)))
            if best is None or rms < best[0]:
                best = (rms, sign, k)
    rms, sign, k = best
    print(f"   planar fit: theta sign {sign:+.0f}, signal order {orders[k]:+.0f} (index {k}), rms vs referee {rms:.2e}", flush=True)
    if rms > 0.05:
        raise RuntimeError(f"no (sign, order) fits the planar referee (best rms {rms:.2e})")
    return sign, k, orders


HEAD = ["set", "case", "phi_deg", "theta_deg", "lambda_um", "order", "theta_lum_deg", "phi_lum_deg",
        "eta_s_lum", "eta_p_lum", "power_s_lum", "power_p_lum", "harmonics", "layers_per_period", "theta_sign"]
HEAD += [f"{k}_{ri}" for k in AMP_KEYS for ri in ("re", "im")]


def sweep_case(fdtd, case_list, case_dict, name, harmonics, layers, bragg_only=False):
    sign, k, orders = planar_sign_fit(fdtd, case_list, case_dict, name, harmonics, layers)
    c0 = case_list[0]
    tau = c0["tau"]
    rows = []
    for c in case_list:
        thetas = np.array([c["theta_bragg_deg"]]) if bragg_only else c["thetas_deg"]
        pairs = lum_pairs(thetas, c["phi_deg"], sign)
        configure_conical(fdtd, case_dict, name, pairs, c["lam_um"] * 1e-6, harmonics, layers)
        res = sept.solve(fdtd, case_dict, False)
        orders_now = np.asarray(res[P["res_orders"]]).ravel()
        hit = np.flatnonzero(orders_now == orders[k])
        if hit.size != 1:
            raise RuntimeError(f"signal order {orders[k]:+.0f} not in this solve's order axis {orders_now.tolist()}")
        kk = int(hit[0])
        r_s, t_s = read_powers(res, "s")
        r_p, t_p = read_powers(res, "p")
        er_s, et_s = sept.total_energy(fdtd, True, "TE")
        er_p, et_p = sept.total_energy(fdtd, True, "TM")
        pow_s = np.asarray(er_s).ravel() + np.asarray(et_s).ravel()
        pow_p = np.asarray(er_p).ravel() + np.asarray(et_p).ravel()
        amps, amp_orders = read_amplitudes(fdtd)
        if amps is not None:
            ha = np.flatnonzero(amp_orders == orders[k])
            ka = int(ha[0]) if ha.size == 1 else None
        for j, th in enumerate(thetas):
            es = (r_s if tau == -1 else t_s)[kk, j]
            ep = (r_p if tau == -1 else t_p)[kk, j]
            row = [c["set"], name, float(c["phi_deg"]), float(th), c["lam_um"], int(orders[k]), float(pairs[j, 0]), float(pairs[j, 1]),
                   float(es), float(ep), float(pow_s[j]), float(pow_p[j]), harmonics, layers, int(sign)]
            if amps is not None and ka is not None:
                for key in AMP_KEYS:
                    v = amps[key][ka, j]
                    row += [float(np.real(v)), float(np.imag(v))]
            else:
                row += [""] * (2 * len(AMP_KEYS))
            rows.append(row)
        # preview against the referee
        ref_s, ref_p = [], []
        for th in thetas:
            R = referee_solve(c["k0"], c["n0"], c["n1"], c["K"], c["Psi"], th * deg, c["phi_deg"] * deg, c["d"], M=10)
            J = R["J_T"] if tau == +1 else R["J_R"]
            ref_s.append(eta_of(J, [1, 0]))
            ref_p.append(eta_of(J, [0, 1]))
        got_s = np.array([r[8] for r in rows[-len(thetas):]])
        got_p = np.array([r[9] for r in rows[-len(thetas):]])
        print(f"   phi {c['phi_deg']:6.1f}: s rms vs referee {np.sqrt(np.mean((got_s - np.array(ref_s)) ** 2)):.2e}, "
              f"p rms {np.sqrt(np.mean((got_p - np.array(ref_p)) ** 2)):.2e}; power sums s {pow_s.min():.6f}..{pow_s.max():.6f}, "
              f"p {pow_p.min():.6f}..{pow_p.max():.6f}   (gate: rms <= 1e-3, sums 1 +- 1e-3)", flush=True)
    return rows


def main(argv):
    harmonics = sept.flag(argv, "--harmonics", None, int)
    layers = sept.flag(argv, "--layers", None, int)
    names = [a for a in argv if not a.startswith("--") and not a.replace(".", "").replace("-", "").isdigit()]
    lumapi = sept.load_lumapi()
    if "--probe" in argv:
        fdtd = sept.open_model(lumapi, sept.FSP)
        print("Lumerical version:", fdtd.version(), flush=True)
        sept.probe(fdtd, "lumerical_R")
        fdtd.close()
        return
    name = names[0]
    case_list = all_cases()[name]
    case_dict, name = model_inputs(case_list)
    fdtd = sept.open_model(lumapi, sept.FSP_UNSLANTED if (name == "unslanted_T" and os.path.isfile(sept.FSP_UNSLANTED)) else sept.FSP)
    print("Lumerical version:", fdtd.version(), flush=True)
    rows, failed = [], False
    try:
        rows = sweep_case(fdtd, case_list, case_dict, name, harmonics, layers, bragg_only="--bragg" in argv)
        if "--bragg" not in argv:
            os.makedirs(os.path.join(ROOT, "lumerical"), exist_ok=True)
            fdtd.save(os.path.join(ROOT, "lumerical", f"rcwa_conical_{cc.file_stem(name)}.fsp"))
    except Exception:
        import traceback
        traceback.print_exc()
        failed = True
    finally:
        if "--bragg" not in argv:
            os.makedirs(RES, exist_ok=True)
            bad = failed or not rows
            out = os.path.join(RES, f"lumerical_conical_{cc.file_stem(name)}{'.FAILED' if bad else ''}.csv")
            with open(out, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(HEAD)
                w.writerows(rows)
            print(f"{len(rows)} rows -> {out}", flush=True)
        else:
            for r in rows:
                print("   ", r[:12])
        try:
            fdtd.close()
        except Exception as e:
            print("close:", e)


if __name__ == "__main__":
    main(sys.argv[1:])
