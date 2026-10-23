"""fpsm_conical -- the two conical F-PSM closed forms of David Brotherton-Ratcliffe's report
"Conical F-PSM closed forms: comparison with vector rigorous coupled-wave analysis and validated
diffraction-efficiency bounds" (Conical_F_PSM_RCWA_Report.pdf, 20 Oct 2026).  Transcribed
2026-10-23; equation numbers are that report's.  Nothing here is derived independently: this is
a transcription, gated by test_conical.py against the report's own Table 2 maxima and Table 3-4
certified bounds, with rcwa_conical.py as the referee.

Three entry points:

  leading(...)         the report's leading model L (Table 2 rows "L"): the four principal
                       states (0,+,s), (0,+,p), (1,tau,s), (1,tau,p) coupled by the first-order
                       vertex only, H = Lambda_P + mu c^(1), identity faces.  Its errors reproduce
                       the printed L rows to four figures (VERIFIED 2026-10-23, sweep_david.py).
  leading_scalar(...)  the two scalar channels of Sec. 4.A, Eq. (58)-(62): g_sigma, g_pi = C_tau
                       g_sigma and the elementary efficiencies (60)-(61) in the fringe
                       (sigma_f, pi_f) basis of Eq. (15).  Equal to leading() at exact Bragg
                       matching; off Bragg the order-1 wave leaves the plane of K and k_0 and the
                       four-state model keeps a small residual cross-coupling that the scalar
                       channels drop, which is why the printed L rows are leading(), not this.
  corrected(...)       the cubic-bulk / quadratic-face model of Appendix B, Eq. (90)-(116): the
                       four principal states, twelve driven companions, the canonical reduction,
                       the face maps and the order-1 Jones matrix J_T or J_R.  U = exp(-i H d) is
                       scipy's matrix exponential, which the report says its validated calculation
                       also uses (Appendix B.C); the quartic radicals of Eq. (103)-(107) are the
                       same function and are not reproduced here.

Units: one length unit throughout (the report uses um); angles in radians.  Conventions and the
exterior basis vectors: basis.py (shared with the referee so the Jones matrices compare in one
basis).  mu = n1 / n0.  Jones matrices are flux-normalised, [output pol, input pol], exterior
(s, p) bases of orders 1 and 0: eta(b) = ||J b||^2.
"""
import numpy as np
from scipy.linalg import expm

from basis import geometry, e_vectors

__all__ = ["leading", "leading_scalar", "corrected", "fringe_rotation", "principal_open"]

# Eq. (94)/(96): Fourier weights f_l^(r) of c, c^2, c^3 (c = cos), and the vertex weights (A_r, B_r)
_F = {1: {1: 0.5, -1: 0.5},
      2: {0: 0.5, 2: 0.25, -2: 0.25},
      3: {1: 0.375, -1: 0.375, 3: 0.125, -3: 0.125}}
_AB = {1: (2.0, 2.0), 2: (1.0, -3.0), 3: (0.0, 4.0)}

D_TOL = 1e-9     # relative size of a principal-companion gap below which the model is refused


def principal_open(k0, n0, K, Psi, theta, phi):
    """True when orders 0 and 1 both propagate (the closed forms need both)."""
    G = geometry(k0, n0, K, Psi, theta, phi, [0, 1])
    return bool(G["gsq"][0] > 0.0 and G["gsq"][1] > 0.0)


# ------------------------------------------------------------------ states and vertices
def _states(tau):
    """Principal set P (Eq. 90) and companion set Q, in the order used throughout."""
    P = [(0, +1, 0), (0, +1, 1), (1, tau, 0), (1, tau, 1)]
    Q = [(0, -1, 0), (0, -1, 1), (1, -tau, 0), (1, -tau, 1)]
    for m in (-1, 2):
        for s in (+1, -1):
            for rho in (0, 1):
                Q.append((m, s, rho))
    return P, Q


def _setup(k0, n0, n1, K, Psi, theta, phi, tau, phi_g):
    """Everything Appendix B.A defines: states, Bloch exponents lambda_i (Eq. 91), the exterior
    electric vectors (92)-(93) and the elementary vertices c^(r)_ij (94)-(96)."""
    P, Q = _states(tau)
    S = P + Q
    nS = len(S)
    G = geometry(k0, n0, K, Psi, theta, phi, [-1, 0, 1, 2])
    beta, q, Ky = G["beta"], G["q"], G["Ky"]
    gam = {m: G["gam"][m + 1] for m in (-1, 0, 1, 2)}
    gsq = {m: G["gsq"][m + 1] for m in (-1, 0, 1, 2)}
    am = {m: G["am"][m + 1] for m in (-1, 0, 1, 2)}
    if not (gsq[0] > 0.0 and gsq[1] > 0.0):
        raise ValueError("a principal order is closed")
    lam = np.array([m * Ky + s * gam[m] for (m, s, _) in S], dtype=complex)          # (91)
    rm = np.array([np.sqrt(gam[m]) for (m, _, _) in S], dtype=complex)
    sig = np.array([s for (_, s, _) in S], dtype=float)
    ev = []
    for (m, s, rho) in S:
        es, ep = e_vectors(am[m], q, gam[m], beta, float(s))                          # (92)-(93)
        ev.append(es if rho == 0 else ep)
    ev = np.array(ev)
    c = {r: np.zeros((nS, nS), dtype=complex) for r in (1, 2, 3)}
    for i, (m, _, _) in enumerate(S):
        for j, (n, _, _) in enumerate(S):
            l = m - n
            T = ev[i, 0] * ev[j, 0] + ev[i, 2] * ev[j, 2]                              # (95)
            V = ev[i, 1] * ev[j, 1]
            pref = sig[i] * beta ** 2 / (2.0 * rm[i] * rm[j]) * np.exp(1j * l * phi_g)
            for r in (1, 2, 3):
                f = _F[r].get(l, 0.0)
                if f:
                    A, B = _AB[r]
                    c[r][i, j] = pref * f * (A * T + B * V)
    return dict(P=P, Q=Q, S=S, lam=lam, sig=sig, c1=c[1], c2=c[2], c3=c[3], Ky=Ky, beta=beta)


def _jones_from_faces(F0, Fd, Ky, d, tau):
    """Eq. (113): J_T = e^{i K_y d} O_T B_T^{-1} E_0,  J_R = O_R B_R^{-1} E_0."""
    E0 = np.zeros((4, 2), dtype=complex)
    E0[0, 0] = E0[1, 1] = 1.0
    if tau == +1:
        return np.exp(1j * Ky * d) * Fd[2:4, :] @ np.linalg.solve(F0, E0)
    BR = np.vstack([F0[0:2, :], Fd[2:4, :]])
    return F0[2:4, :] @ np.linalg.solve(BR, E0)


# ------------------------------------------------------------------ leading model L
def leading(k0, n0, n1, K, Psi, theta, phi, d, tau, phi_g=0.0):
    """The report's leading model: H = Lambda_P + mu c^(1)_PP, U = exp(-i H d), faces F_0 = I,
    F_d = U.  Returns the order-1 Jones matrix (2 x 2)."""
    mu = n1 / n0
    s = _setup(k0, n0, n1, K, Psi, theta, phi, tau, phi_g)
    H = np.diag(s["lam"][:4]) + mu * s["c1"][:4, :4]
    U = expm(-1j * H * d)
    return _jones_from_faces(np.eye(4, dtype=complex), U, s["Ky"], d, tau)


# ------------------------------------------------------------------ scalar channels, Eq. (58)-(62)
def fringe_rotation(k0, n0, K, Psi, theta, phi, m, sigma):
    """R[a, rho] = f_a . e_rho for order m travelling in direction sigma: the real orthogonal map
    from the exterior (s, p) basis to the fringe (sigma_f, pi_f) basis of Eq. (15)-(16).
    s_f = e_xi x t with t the unit fringe-parallel part of k, so s_f is normal to the plane of K
    and k; at Bragg matching it is the same vector for orders 0 and 1."""
    G = geometry(k0, n0, K, Psi, theta, phi, [m])
    beta, q, am, gam = G["beta"], G["q"], G["am"][0], G["gam"][0]
    if G["gsq"][0] <= 0.0:
        raise ValueError("fringe basis asked for a closed order")
    k = np.array([am, sigma * gam.real, q])
    et = np.array([np.cos(Psi), np.sin(Psi), 0.0])
    exi = np.array([-np.sin(Psi), np.cos(Psi), 0.0])
    p = float(k @ et)
    rho = np.hypot(p, q)
    if rho < 1e-12 * beta:
        raise ValueError("fringe basis singular: k parallel to the grating vector")
    t = (p * et + q * np.array([0.0, 0.0, 1.0])) / rho
    sf = np.cross(exi, t)
    pf = np.cross(k / beta, sf)
    es, ep = e_vectors(am, q, gam, beta, sigma)
    return np.array([[sf @ es.real, sf @ ep.real],
                     [pf @ es.real, pf @ ep.real]])


def leading_scalar(k0, n0, n1, K, Psi, theta, phi, d, tau):
    """Sec. 4.A, Eq. (58)-(62).  Returns dict: g_sigma, g_pi, C_tau, delta, eta_sigma, eta_pi, and
    CtC = R_0^T diag(eta_sigma, eta_pi) R_0, the efficiency operator in the order-0 exterior (s, p)
    basis (output-side rotations are orthogonal and drop out of C^H C)."""
    G = geometry(k0, n0, K, Psi, theta, phi, [0, 1])
    beta, alpha, q, Kx, Ky = G["beta"], G["alpha"], G["q"], G["Kx"], G["Ky"]
    if not (G["gsq"][0] > 0.0 and G["gsq"][1] > 0.0):
        raise ValueError("a principal order is closed")
    g0, g1 = G["gam"][0].real, G["gam"][1].real
    g_sigma = k0 * k0 * n0 * n1 / (2.0 * np.sqrt(g0 * g1))                        # (58)
    C_tau = (alpha * (alpha + Kx) + tau * g0 * g1 + q * q) / beta ** 2              # (59)
    g_pi = C_tau * g_sigma
    delta = 0.5 * (g0 - Ky - tau * g1)                                              # (59)
    out = dict(g_sigma=g_sigma, g_pi=g_pi, C_tau=C_tau, delta=delta)
    for name, g in (("sigma", g_sigma), ("pi", g_pi)):
        if tau == +1:
            Om = np.sqrt(delta ** 2 + g ** 2)
            eta = g * g * np.sin(d * Om) ** 2 / (Om * Om)                           # (60)
        else:
            w = np.sqrt(complex(g * g - delta * delta))                             # (62)
            if w == 0:
                amp = g * d / (1.0 + 1j * delta * d)
            else:
                amp = (g * np.sinh(w * d) / w) / (np.cosh(w * d) + 1j * delta * np.sinh(w * d) / w)   # (61)
            eta = float(abs(amp) ** 2)
        out["eta_" + name] = float(eta)
    R0 = fringe_rotation(k0, n0, K, Psi, theta, phi, 0, +1.0)
    out["CtC"] = R0.T @ np.diag([out["eta_sigma"], out["eta_pi"]]) @ R0
    return out


# ------------------------------------------------------------------ corrected model, Appendix B
def corrected(k0, n0, n1, K, Psi, theta, phi, d, tau, phi_g=0.0, return_parts=False):
    """Appendix B.  Returns the flux-normalised order-1 Jones matrix J (2 x 2): J_T for tau = +1,
    J_R for tau = -1.  Raises when a principal order is closed or a principal-companion gap
    vanishes (the report regularises neither)."""
    mu = n1 / n0
    s = _setup(k0, n0, n1, K, Psi, theta, phi, tau, phi_g)
    lam, sig, c1, c2, c3, Ky, beta = s["lam"], s["sig"], s["c1"], s["c2"], s["c3"], s["Ky"], s["beta"]
    nP, nS = 4, len(s["S"])
    iP, iQ = slice(0, nP), slice(nP, nS)

    # canonical bulk reduction, Eq. (97)-(102)
    lamP, lamQ = lam[iP], lam[iQ]
    dmat = lamP[:, None] - lamQ[None, :]                                              # d_iv, [i in P, v in Q]
    if np.min(np.abs(dmat)) < D_TOL * beta:
        raise ValueError("a principal-companion gap vanishes: enlarge the core")
    X = c1[iQ, iP] / dmat.T                                                           # X_vi, [v, i]
    Y = -c1[iP, iQ] / dmat                                                            # Y_iv, [i, v]
    X2 = (c2[iQ, iP] + c1[iQ, iQ] @ X - X @ c1[iP, iP]) / dmat.T                       # (98)
    Y2 = (-c2[iP, iQ] - c1[iP, iP] @ Y + Y @ c1[iQ, iQ]) / dmat                        # (99)
    inv_d = 1.0 / dmat                                                                # [i, v]
    h2 = c2[iP, iP] + 0.5 * ((c1[iP, iQ] * inv_d) @ c1[iQ, iP] + c1[iP, iQ] @ (c1[iQ, iP] * inv_d.T))   # (100)
    h3 = c3[iP, iP] + 0.5 * (c1[iP, iQ] @ X2 + c2[iP, iQ] @ X - Y2 @ c1[iQ, iP] - Y @ c2[iQ, iP])      # (101)
    H = np.diag(lamP) + mu * c1[iP, iP] + mu ** 2 * h2 + mu ** 3 * h3                 # (102)
    U = expm(-1j * H * d)

    # companion propagators and projectors, Eq. (108)-(109); only the decaying branch of each
    # companion is ever formed (the other would overflow for a closed order at these thicknesses)
    sQ = sig[iQ]
    Pi_p = np.diag((sQ > 0).astype(float)).astype(complex)
    Pi_m = np.diag((sQ < 0).astype(float)).astype(complex)
    E_p = np.diag(np.array([np.exp(-1j * l * d) if sg > 0 else 0.0 for l, sg in zip(lamQ, sQ)], dtype=complex))
    E_m = np.diag(np.array([np.exp(+1j * l * d) if sg < 0 else 0.0 for l, sg in zip(lamQ, sQ)], dtype=complex))
    if tau == +1:
        # the backward quartet B = Q[0:4] propagates jointly through exp(+i H_b d), Eq. (109)
        iB = [nP + k for k in range(4)]
        rest = [k for k in range(nS) if k not in iB]
        Hb = np.diag(lam[iB]) + mu * c1[np.ix_(iB, iB)] + mu ** 2 * c2[np.ix_(iB, iB)]
        for a, v in enumerate(iB):
            for b, w in enumerate(iB):
                ssum = 0.0
                for r in rest:
                    gv, gw = lam[v] - lam[r], lam[w] - lam[r]
                    if min(abs(gv), abs(gw)) < D_TOL * beta:
                        # e.g. normal incidence on a transmission grating at exact Bragg, where the
                        # (-1,-) companion is degenerate with (0,-); not regularised, by the report's rule
                        raise ValueError("a companion-companion gap in H_b vanishes: enlarge the core")
                    ssum += c1[v, r] * c1[r, w] * (1.0 / gv + 1.0 / gw)
                Hb[a, b] += mu ** 2 * 0.5 * ssum
        E_m[0:4, 0:4] = expm(+1j * Hb * d)

    # faces, Eq. (110)-(112)
    Af = np.eye(nP, dtype=complex) + 0.5 * mu ** 2 * (Y @ X)
    F0 = Af - mu ** 2 * Y @ (Pi_p @ X + E_m @ X @ U)
    Fd = Af @ U - mu ** 2 * Y @ (E_p @ X + Pi_m @ X @ U)
    J = _jones_from_faces(F0, Fd, Ky, d, tau)
    if return_parts:
        return J, dict(H=H, U=U, X=X, Y=Y, F0=F0, Fd=Fd, c1=c1, c2=c2, c3=c3, lam=lam, S=s["S"])
    return J
