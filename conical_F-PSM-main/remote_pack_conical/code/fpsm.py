"""fpsm -- the two F-PSM closed-form efficiency formulae (lossless, TE/sigma).

Source: "The two simple F-PSM efficiency formulae -- an implementation note for
lossless TE (sigma) planar gratings" (``../FPSM_Simple_Implementation.pdf``).  Equation numbers below are
that note's.  Nothing here is derived independently: this is a transcription,
checked against the note's Section 6 test vectors by ``test_fpsm_vectors.py``.

Conventions (note, "Units and conventions"): exp(+i w t); every length in ONE
unit (the formulae are scale-free); angles passed in RADIANS to this module;
K is SIGNED (eq. 2) and must never be passed through abs().

    theta_c = theta_inc + psi        angle between the input beam and K
    K       = 4 pi n0 cos(theta_r) / lam_r                        (eq. 2)

Public entry point: ``efficiency`` (vectorised over any broadcastable inputs).
It never refuses to return a number outside the stated working domain (note
Section 5) -- it reports the domain flags alongside.  It DOES refuse where the
note says "hard stop": a closed principal channel returns NaN, never a number.
"""
from __future__ import annotations

import numpy as np

__all__ = ["efficiency", "fpsm_detail", "K_from_recording", "recover_K_sign"]


def K_from_recording(n0, lam_r, theta_r):
    """Signed grating-vector magnitude, eq. (2).  theta_r in radians."""
    return 4.0 * np.pi * n0 * np.cos(theta_r) / lam_r


def _outgoing_root(gsq):
    """Gamma of eq. (38): sqrt(gsq) if open, -i sqrt(-gsq) if closed.
    Written as the explicit two-branch form the note asks for."""
    gsq = np.asarray(gsq, dtype=float)
    mag = np.sqrt(np.abs(gsq))
    return np.where(gsq > 0.0, mag + 0j, -1j * mag)


def _geometry(k0, n0, K, psi, theta_c):
    beta = k0 * n0                                            # eq. (3)
    Kx, Ky = K * np.sin(psi), K * np.cos(psi)                 # eq. (4)
    alpha = beta * np.sin(theta_c - psi)                      # eq. (5)
    gsq = {m: beta ** 2 - (alpha + m * Kx) ** 2 for m in (-1, 0, 1, 2)}   # eq. (6)
    return beta, Kx, Ky, gsq


def recover_K_sign(k0, n0, K_abs, psi, theta_c):
    """Note Section 2.2: when only the fringe spacing is known, evaluate both
    signs and both residuals and keep the smallest of the four.  Scalar inputs.
    Returns the signed K."""
    best = None
    for s in (+1.0, -1.0):
        _, _, Ky, gsq = _geometry(k0, n0, s * K_abs, psi, theta_c)
        if gsq[0] <= 0 or gsq[1] <= 0:
            continue
        g0, g1 = np.sqrt(gsq[0]), np.sqrt(gsq[1])
        r = min(abs(g0 - g1 - Ky), abs(g0 + g1 - Ky))
        if best is None or r < best[0]:
            best = (r, s * K_abs)
    if best is None:
        raise ValueError("principal channel closed for both signs of K")
    return best[1]


TM_ABSENT = ("F-PSM exists for lossless TE only: the implementation note of 16 September 2026 and the "
             "Journal of Optics submission of 6 October 2026 both restrict it so, and the paper states "
             "that complex coefficients in the lossless formulae are NOT an absorption theory. No TM or "
             "mixed-grating sheet had been received as of 2026-10-16. Nothing is substituted.")


def _require_te_lossless(pol, chi0, chi1):
    """Fail loudly (never silently fall back to TE) when a polarisation or loss F-PSM does not have."""
    if pol != "TE":
        raise NotImplementedError(f"F-PSM pol={pol!r}: " + TM_ABSENT)
    if np.any(np.asarray(chi0, dtype=float) != 0.0) or np.any(np.asarray(chi1, dtype=float) != 0.0):
        raise NotImplementedError("F-PSM with absorption (chi0 or chi1 != 0): " + TM_ABSENT)


def fpsm_detail(k0, n0, n1, K, psi, theta_c, d, pol="TE", chi0=0.0, chi1=0.0):
    """Full evaluation.  All inputs broadcast; angles in radians.  ``pol``, ``chi0``, ``chi1`` are
    accepted so every engine shares one signature; anything but TE lossless raises (2026-10-16).

    Returns a dict of arrays: ``eta`` (NaN where the principal channel is
    closed), ``case`` (+1 transmission, -1 reflection, 0 not applicable),
    ``in_domain`` (all six Section-5 conditions hold), ``domain_fail`` (bitmask:
    bit j-1 set when condition Tj / Rj fails), and the intermediates the note's
    Section 6 tabulates (g, ge, delta, gamma0, gamma1, Hs, Bm, Bp).
    """
    _require_te_lossless(pol, chi0, chi1)
    k0, n0, n1, K, psi, theta_c, d = np.broadcast_arrays(
        *[np.asarray(v, dtype=float) for v in (k0, n0, n1, K, psi, theta_c, d)])
    beta, Kx, Ky, gsq = _geometry(k0, n0, K, psi, theta_c)

    applicable = (gsq[0] > 0.0) & (gsq[1] > 0.0)              # hard stop, eq. (7)
    g0 = np.sqrt(np.where(applicable, gsq[0], 1.0))
    g1 = np.sqrt(np.where(applicable, gsq[1], 1.0))
    gm1sq, gp2sq = gsq[-1], gsq[2]                            # kept as SQUARES

    V = k0 ** 2 * n0 * n1                                     # eq. (9)
    V0 = 0.5 * k0 ** 2 * n1 ** 2
    V2 = 0.25 * k0 ** 2 * n1 ** 2

    RT = np.abs(g0 - g1 - Ky)                                 # eq. (8)
    RR = np.abs(g0 + g1 - Ky)
    is_T = RT < RR

    with np.errstate(divide="ignore", invalid="ignore"):
        T = _transmission(g0, g1, Ky, gm1sq, gp2sq, V, V0, V2, d, n0, n1, beta)
        R = _reflection(g0, g1, Ky, gm1sq, gp2sq, V, V0, V2, d, n0, n1, beta)

    out = {}
    for key in T:
        out[key] = np.where(is_T, T[key], R[key])
    out["Hs"] = np.where(is_T, np.nan, R["Hs"])
    out["Bm"] = np.where(is_T, np.nan, R["Bm"])
    out["Bp"] = np.where(is_T, np.nan, R["Bp"])
    out["eta"] = np.where(applicable, out["eta"], np.nan)
    out["case"] = np.where(applicable, np.where(is_T, 1, -1), 0)
    out["gamma0"], out["gamma1"] = g0, g1
    out["in_domain"] = applicable & (out["domain_fail"] == 0)
    return out


def _transmission(g0, g1, Ky, gm1sq, gp2sq, V, V0, V2, d, n0, n1, beta):
    g = V / (2.0 * np.sqrt(g0 * g1))                          # eq. (10)
    L = g0 + g1 - Ky                                          # eq. (11)
    M = g0 + g1 + Ky
    r0, r1, r01 = g ** 2 / L ** 2, g ** 2 / M ** 2, g ** 2 / (L * M)      # eq. (12)

    xm0, xm1 = g0 + Ky, g1 + 2.0 * Ky                         # eq. (13)
    xp0, xp1 = g0 - 2.0 * Ky, g1 - Ky                         # eq. (14)
    Dm0, Dm1 = xm0 ** 2 - gm1sq, xm1 ** 2 - gm1sq
    Dp0, Dp1 = xp0 ** 2 - gp2sq, xp1 ** 2 - gp2sq

    p0 = g0 + V0 / (2 * g0) - g ** 2 / L + V ** 2 / (2 * g0 * Dm0)        # eq. (15)
    p1 = Ky + g1 + V0 / (2 * g1) - g ** 2 / M + V ** 2 / (2 * g1 * Dp1)   # eq. (16)
    p, delta = 0.5 * (p0 + p1), 0.5 * (p0 - p1)               # eq. (17)

    eb = g ** 2 / 4 * (1 / (g1 * L) + 1 / (g0 * M)) + g ** 2 / (L * M)    # eq. (18)
    em = -V0 / 4 * ((1 / g0) * (1 / M + 1 / (2 * g0))
                    + (1 / g1) * (1 / L + 1 / (2 * g1)))                  # eq. (19)
    eq = (-V ** 2 / 4 * ((xm0 + xm1) / (g0 * Dm0 * Dm1)
                         + (xp0 + xp1) / (g1 * Dp0 * Dp1))
          + V2 / 2 * (1 / Dm0 + 1 / Dm1 + 1 / Dp0 + 1 / Dp1)
          - V ** 2 / (2 * g0 * M * Dm0)
          - V ** 2 / (2 * g1 * L * Dp1))                                  # eq. (20)
    ge = g * (1 + eb + em + eq)

    pbar = 0.5 * (Ky - g0 - g1)                               # eq. (21)
    dbar = 0.5 * (-g0 + g1 - Ky)

    Om = np.sqrt(delta ** 2 + ge ** 2)                        # eq. (22)
    Omb = np.sqrt(dbar ** 2 + g ** 2)                         # eq. (23)
    F = np.where(Om == 0, d, np.sin(Om * d) / np.where(Om == 0, 1.0, Om))
    Fb = np.where(Omb == 0, d, np.sin(Omb * d) / np.where(Omb == 0, 1.0, Omb))
    C, Cb = np.cos(Om * d), np.cos(Omb * d)
    A = ge * F

    J = ((r0 + r1) * (Cb * C - dbar * delta * Fb * F)
         + 1j * (r1 - r0) * (Cb * delta * F + dbar * Fb * C))             # eq. (24)
    amp = ((1 - 0.5 * (r0 + r1)) * A
           + np.exp(1j * (pbar - p) * d) * (A * J + r01 * g * Fb * (1 - 2 * A ** 2)))  # eq. (25)
    eta = np.abs(amp) ** 2                                    # eq. (26)

    # Section 5.1 working domain, T1..T6
    mLM = np.minimum(np.abs(L), np.abs(M))
    closed = (gm1sq < 0) & (gp2sq < 0)
    decay = np.minimum(np.sqrt(np.abs(gm1sq)), np.sqrt(np.abs(gp2sq)))
    conds = [np.abs(n1) / n0 <= 0.03,
             np.minimum(g0, g1) / beta >= 0.25,
             np.abs(g) / mLM <= 0.1,
             np.abs(g0 - Ky - g1) / (2 * mLM) <= 0.1,
             closed,
             closed & (d * decay >= 8.0)]
    fail = sum((~c).astype(int) << j for j, c in enumerate(conds))
    return dict(eta=eta, g=g, ge=ge, delta=delta, domain_fail=fail)


def _reflection(a, b, Ky, gm1sq, gp2sq, V, V0, V2, d, n0, n1, beta):
    g = V / (2.0 * np.sqrt(a * b))                            # eq. (28)
    L = a - b - Ky                                            # eq. (29)
    M = a - b + Ky

    xm0, xm1 = a + Ky, -b + 2.0 * Ky                          # eq. (30)
    xp0, xp1 = a - 2.0 * Ky, -b - Ky                          # eq. (31)
    Dm0, Dm1 = xm0 ** 2 - gm1sq, xm1 ** 2 - gm1sq
    Dp0, Dp1 = xp0 ** 2 - gp2sq, xp1 ** 2 - gp2sq

    p0 = a + V0 / (2 * a) + g ** 2 / L + V ** 2 / (2 * a * Dm0)           # eq. (32)
    p1 = Ky - b - V0 / (2 * b) + g ** 2 / M - V ** 2 / (2 * b * Dp1)      # eq. (33)
    p, delta = 0.5 * (p0 + p1), 0.5 * (p0 - p1)               # eq. (34)

    eb = g ** 2 / 4 * (1 / (b * L) - 1 / (a * M)) - g ** 2 / (L * M)      # eq. (35)
    em = -V0 / 4 * ((1 / a) * (1 / M + 1 / (2 * a))
                    - (1 / b) * (1 / L - 1 / (2 * b)))                    # eq. (36)
    eq = (-V ** 2 / 4 * ((xm0 + xm1) / (a * Dm0 * Dm1)
                         - (xp0 + xp1) / (b * Dp0 * Dp1))
          + V2 / 2 * (1 / Dm0 + 1 / Dm1 + 1 / Dp0 + 1 / Dp1)
          - V ** 2 / (2 * a * M * Dm0)
          + V ** 2 / (2 * b * L * Dp1))                                   # eq. (37)
    ge = g * (1 + eb + em + eq)

    Gm, Gp = _outgoing_root(gm1sq), _outgoing_root(gp2sq)     # eq. (38)
    Hs = (g ** 2 / 2 * (1 / L ** 2 + 1 / M ** 2)
          + V ** 2 * (xm0 ** 2 + Gm ** 2) / (4 * a * Gm * Dm0 ** 2)
          + V ** 2 * (xp1 ** 2 + Gp ** 2) / (4 * b * Gp * Dp1 ** 2))      # eq. (39)
    Bm = V ** 2 * np.exp(-1j * (p + Ky + Gm) * d) / (4 * a * Gm * (xm0 + Gm) ** 2)   # eq. (40)
    Bp = V ** 2 * np.exp(-1j * (2 * Ky + Gp - p) * d) / (4 * b * Gp * (xp1 - Gp) ** 2)

    q = np.sqrt(ge ** 2 - delta ** 2 + 0j)                    # eq. (41), complex, Re(q) >= 0
    # eq. (42) divided through by cosh(qd): identical algebra, no overflow for
    # thick plates.  tq = tanh(qd)/q -> d as q -> 0; sech(qd) -> 0 as Re(qd) -> inf.
    qd = q * d
    e1, e2 = np.exp(-qd), np.exp(-2.0 * qd)
    sech = 2.0 * e1 / (1.0 + e2)
    tq = np.where(q == 0, d + 0j, (1.0 - e2) / (1.0 + e2) / np.where(q == 0, 1.0, q))
    eta = np.abs(ge * tq / ((1 + Hs) * (1 + 1j * delta * tq) - (Bm + Bp) * sech)) ** 2

    # Section 5.2 working domain, R1..R6
    G = np.minimum.reduce([np.abs(L), np.abs(M),
                           np.abs(a + Ky + Gm), np.abs(a + Ky - Gm),
                           np.abs(-b - Ky + Gp), np.abs(-b - Ky - Gp)])   # eq. (44)
    gmax = np.abs(V) / (2 * np.sqrt(np.minimum(a, b)
                                    * np.minimum.reduce([a, b, np.abs(Gm), np.abs(Gp)])))  # eq. (45)
    conds = [np.abs(n1) / n0 <= 0.03,
             np.minimum(a, b) / beta >= 0.25,
             np.minimum(np.abs(Gm), np.abs(Gp)) / beta >= 0.2,
             G / beta >= 0.25,
             gmax / G <= 0.1,
             np.abs(a + b - Ky) / (2 * G) <= 0.1]
    fail = sum((~c).astype(int) << j for j, c in enumerate(conds))
    return dict(eta=eta, g=g, ge=ge, delta=delta, domain_fail=fail, Hs=Hs, Bm=Bm, Bp=Bp)


def efficiency(k0, n0, n1, K, psi, theta_c, d, pol="TE", chi0=0.0, chi1=0.0):
    """eta in [0, 1]; NaN where the principal channel is closed ("not
    applicable", note eq. 7).  Angles in radians, lengths in one unit.
    pol != "TE" or any absorption raises NotImplementedError (see TM_ABSENT)."""
    return fpsm_detail(k0, n0, n1, K, psi, theta_c, d, pol=pol, chi0=chi0, chi1=chi1)["eta"]
