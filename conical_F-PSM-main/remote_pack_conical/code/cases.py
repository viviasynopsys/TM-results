"""The benchmark gratings.  One definition feeds F-PSM, Kogelnik, RCWA and the
OpticStudio bench, so the four can never drift apart.

A grating is given the way the Ansys articles give it: two collimated recording
beams INSIDE the medium, as unit vectors (t, n) = (tangential, normal) with the
normal pointing from the entry face into the plate.  Grating vector
K_vec = beta_r (k1 - k2); replay near k1 diffracts into k2.

F-PSM geometry (note Section 1.3): K_vec = K (-sin psi, cos psi).  We take
K = +|K_vec| and psi = atan2(-K_t, K_n).  (K, psi) -> (-K, psi + 180 deg) leaves
K_x, K_y and theta_inc untouched, so this is the same grating as the note's
signed-K form; ``test_engines.py`` confirms the note's T-a/T-b vectors survive it.

Replay angle theta (the articles' axis): input k = beta (sin theta, cos theta),
theta_c = theta + psi.  OpticStudio bench: plane of incidence y-z, (t, n) ->
(y, z); a Source Ray tilted tau about X points along (0, -sin tau, cos tau), so
tau = -theta.  Holo type 2: beam 1 DIVERGES from P1 = -k1 * far, beam 2
CONVERGES to P2 = +k2 * far.
"""
import numpy as np

deg = np.pi / 180.0
FAR = 1.0e8   # mm; OpticStudio treats a far construction point as a collimated beam


def _unit(angle_deg, forward=True):
    a = angle_deg * deg
    s = 1.0 if forward else -1.0
    return np.array([np.sin(a), s * np.cos(a)])


def _case(title, source, n0, n1, d_um, lam_um, k1, k2, theta_deg, lam_sweep_um, lam_r_um=None):
    lam_r_um = lam_um if lam_r_um is None else lam_r_um
    Kvec = (2 * np.pi * n0 / (lam_r_um * 1e3)) * (k1 - k2)          # nm^-1
    K = float(np.hypot(*Kvec))
    psi = float(np.arctan2(-Kvec[0], Kvec[1]))
    theta_bragg = float(np.arctan2(k1[0], k1[1]))
    return dict(
        title=title, source=source, n0=n0, n1=n1, d_um=d_um, d_nm=d_um * 1e3,
        lam_um=lam_um, lam_nm=lam_um * 1e3, lam_r_um=lam_r_um,
        K=K, psi=psi, theta_bragg_deg=theta_bragg / deg,
        theta_deg=np.asarray(theta_deg, dtype=float),
        lam_sweep_nm=np.asarray(lam_sweep_um, dtype=float) * 1e3,
        P1=[0.0, float(-k1[0] * FAR), float(-k1[1] * FAR)],
        P2=[0.0, float(k2[0] * FAR), float(k2[1] * FAR)],
        zemax_tilts_deg=[-float(t) for t in theta_deg],
        zemax_tilt_bragg_deg=-theta_bragg / deg,
        zemax_lams_um=[float(v) for v in lam_sweep_um],
    )


CASES = {
    # Ansys Lumerical article 37118679173267 "Volume Holographic Grating": n0 = 1.5,
    # dn = 0.04, 6 um, 550 nm, beams at 0 and 28 deg in the medium from opposite
    # sides, background index 1.5, S-pol; its figure plots Rs order +1, theta -50..30.
    "lumerical_R": _case(
        "Reflection VHG of the Ansys Lumerical article", "optics.ansys.com 37118679173267",
        1.5, 0.04, 6.0, 0.55, _unit(0.0), _unit(28.0, forward=False),
        np.arange(-50.0, 30.0 + 1e-9, 0.5), np.arange(0.45, 0.65 + 1e-9, 0.002)),
    # Same article, "T" variant: sign of k2z flipped.  Neighbour order m = 2 is OPEN
    # (condition T5 fails) -- deliberately outside F-PSM's stated domain.
    "lumerical_T": _case(
        "Transmission variant of the Ansys Lumerical article", "optics.ansys.com 37118679173267",
        1.5, 0.04, 6.0, 0.55, _unit(0.0), _unit(28.0),
        np.arange(-20.0, 50.0 + 1e-9, 0.5), np.arange(0.40, 0.70 + 1e-9, 0.003)),
    # Ansys Zemax article 42661816933139, attachment test_example1.zmx (RCWA DLL demo):
    # n0 = 1.5, dn = 0.05, 5 um, 0.55 um, P1 = (0, 1, -1), P2 = (0, 0, -1), n1 = n2 = 1.5,
    # source ray tilted 45 deg.  dn/n0 = 0.033 is just outside F-PSM condition R1.
    "zemax_kba_R": _case(
        "Reflection VHG of the Ansys Zemax RCWA article (test_example1)", "optics.ansys.com 42661816933139",
        1.5, 0.05, 5.0, 0.55, _unit(-45.0), _unit(0.0, forward=False),
        np.arange(-70.0, 25.0 + 1e-9, 0.5), np.arange(0.45, 0.65 + 1e-9, 0.002)),
    # the F-PSM note, Section 6.1 case T-a: unslanted transmission, inside the stated domain.
    "unslanted_T": _case(
        "Unslanted transmission grating, the note's test vector T-a", "FPSM_Simple_Implementation.pdf 6.1",
        1.5, 0.03, 10.0, 0.50, _unit(33.8627180249), _unit(-33.8627180249),
        np.arange(24.0, 44.0 + 1e-9, 0.1), np.arange(0.40, 0.60 + 1e-9, 0.002)),
    # ADDED 2026-09-18 for the OpticStudio RCWA DLL only.  The DLL slices one z-PERIOD of the
    # grating into "Layer per period" layers; an exactly unslanted grating (K_z = 0, infinite
    # z-period) kills the OpticStudio server, and rotations so small that the 10 um slab is thinner
    # than one z-period (0.01 / 0.1 / 1 deg) return 0.88 / 0.49 / 0.05 at Bragg against 0.586
    # (results/zemax_rcwa_probe5.txt).  This is the note's T-a grating with BOTH recording beams
    # turned by 5 deg (same |K|, n1, d; psi = -95 deg; slab = 2.9 z-periods): the nearest grating
    # the DLL can represent.  It is scored against the referee on this same rotated geometry; it is
    # NOT the note's vector T-a and is labelled as a substitute wherever it appears.
    "unslanted_T_rot5": _case(
        "T-a rotated 5° (DLL substitute for the unslanted grating)", "execution.md step 3 deviation, 2026-09-18",
        1.5, 0.03, 10.0, 0.50, _unit(33.8627180249 + 5.0), _unit(-33.8627180249 + 5.0),
        np.arange(29.0, 49.0 + 1e-9, 0.1), np.arange(0.40, 0.60 + 1e-9, 0.002)),
}
