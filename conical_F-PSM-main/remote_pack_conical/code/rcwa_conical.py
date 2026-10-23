"""rcwa_conical -- vector rigorous coupled-wave referee for a slanted sinusoidal volume grating
at conical incidence (any polar angle theta, any azimuth phi).

Written 2026-10-23 to test David Brotherton-Ratcliffe's conical F-PSM closed forms
(Conical_F_PSM_RCWA_Report.pdf, 20 Oct 2026; "the report").  The 4N first-order system is derived
here from Maxwell's curl equations in the report's F-PSM coordinates (Sec. 2.A), not transcribed
from the report's Eq. (10)-(11) nor from Moharam & Gaylord (1983).  test_conical.py ties it to the
planar two-polarisation referee of 2026-10-16 (rcwa.py) at phi = 0, to power conservation at
conical points, and to the Moharam-Gaylord Fig. 5 geometry the report recomputes.

Conventions: basis.py.  With fields ~ exp(-i(alpha_m x' + q z)) exp(i m K_y y) F_m(y) and
g = Z0 H, the curl equations give, for the BULK tangential amplitudes F = (E_x, E_z, g_x, g_z)
(every component a vector over the N = 2M+1 retained orders),

    i F' = M F,
    g_y = (q E_x - alpha E_z)/k0                                  (exact, per order)
    E_y = (1/k0) A (alpha g_z - q g_x)                            (A = Toeplitz of 1/eps)
    i E_x' = L E_x + alpha E_y - k0 g_z
    i E_z' = L E_z + q E_y + k0 g_x
    i g_x' = L g_x + alpha g_y + k0 Eps E_z                       (Eps = Toeplitz of eps)
    i g_z' = L g_z + q g_y - k0 Eps E_x

with L = diag(m K_y), alpha = diag(alpha_m).  Modes F ~ W exp(-i lam y); those with Im(lam) > 0
grow toward +y and are anchored at y = d, the rest at y = 0, so nothing overflows.  Both faces are
matched to the exterior plane-wave set (every retained order, both polarisations, both directions)
in one linear solve with two right-hand sides (s and p incidence).

Note on the report's Eq. (11): as printed, P_11 = q A cal(A)/k0 and P_22 = -q A cal(A)/k0 have the
same operator order, whereas the elimination of E_y gives alpha.A.(...) in the E_x row and
q A.alpha in the E_z row, i.e. P_22 = -q cal(A) A/k0.  The two differ by the commutator
[A, cal(A)] when q != 0.  test_conical.py gate G5 evaluates both orders; the printed one breaks
power conservation, the derived one conserves it.  (A = diag(a_m), cal(A) = Toeplitz(1/eps).)
"""
import numpy as np

from basis import geometry, tangential_fields, material_harmonics

__all__ = ["solve", "jones_error", "eta_of"]


def _toeplitz(coef, N):
    """T[m, n] = coef[m - n] for a dict keyed by signed h; zero where absent."""
    T = np.zeros((N, N), dtype=complex)
    for m in range(N):
        for n in range(N):
            T[m, n] = coef.get(m - n, 0.0)
    return T


def solve(k0, n0, n1, K, Psi, theta, phi, d, M=9, phi_g=0.0, eps_h=None, inv_eps_h=None,
          p22_order="derived"):
    """One conical point.  Scalars; angles in radians; one length unit.

    eps_h / inv_eps_h: optional dicts (signed h -> coefficient) that REPLACE the sinusoidal-index
    material -- used for the report's Moharam-Gaylord check, whose grating is a sinusoidal
    PERMITTIVITY.  n1 is ignored when they are given.

    Returns a dict:
      ms, gam, open   per-order index, normal wavenumber, open flag
      DE_R, DE_T      [order, incident pol] efficiencies (fraction of incident power), s then p
      sum_s, sum_p    total open-order power for s and p incidence (1 when lossless)
      J_T, J_R        order-1 flux-normalised Jones matrices, [output pol, input pol] in basis.py's e-vectors
      t, r            raw electric amplitudes [order, output pol, input pol] at y = d and y = 0
    """
    N = 2 * M + 1
    ms = np.arange(-M, M + 1)
    G = geometry(k0, n0, K, Psi, theta, phi, ms)
    beta, q, Ky, am, gam, gsq = G["beta"], G["q"], G["Ky"], G["am"], G["gam"], G["gsq"]
    if gsq[M] <= 0.0:
        raise ValueError("the incident order is closed: theta is not a propagating direction")

    if eps_h is None:
        eps_h, inv_eps_h = material_harmonics(n0, n1, phi_g, N - 1)
    elif inv_eps_h is None:
        raise ValueError("eps_h given without inv_eps_h")
    Eps = _toeplitz(eps_h, N)
    A = _toeplitz(inv_eps_h, N)

    L = np.diag(ms * Ky).astype(complex)
    Al = np.diag(am).astype(complex)
    I = np.eye(N, dtype=complex)
    Z = np.zeros((N, N), dtype=complex)
    if p22_order == "derived":
        P22 = (q / k0) * (A @ Al)
    elif p22_order == "printed":
        P22 = (q / k0) * (Al @ A)
    else:
        raise ValueError("p22_order must be 'derived' or 'printed'")
    Mx = np.block([
        [L, Z, -(q / k0) * (Al @ A), (Al @ A @ Al) / k0 - k0 * I],
        [Z, L, -(q * q / k0) * A + k0 * I, P22],
        [(q / k0) * Al, -(Al @ Al) / k0 + k0 * Eps, L, Z],
        [(q * q / k0) * I - k0 * Eps, -(q / k0) * Al, Z, L],
    ])
    lam, W = np.linalg.eig(Mx)
    up = lam.imag > 0.0
    s0 = np.ones(4 * N, dtype=complex)
    sd = np.ones(4 * N, dtype=complex)
    s0[up] = np.exp(1j * lam[up] * d)          # exp(-i lam (0 - d)), |.| < 1
    sd[~up] = np.exp(-1j * lam[~up] * d)       # exp(-i lam d), |.| <= 1
    phase_d = np.tile(np.exp(1j * ms * Ky * d), 4)

    Vm = np.zeros((4 * N, 2 * N), dtype=complex)     # exterior waves travelling toward -y (reflected)
    Vp = np.zeros((4 * N, 2 * N), dtype=complex)     # toward +y (transmitted)
    for j in range(N):
        rows = [j, N + j, 2 * N + j, 3 * N + j]
        Vm[rows, 2 * j:2 * j + 2] = tangential_fields(am[j], q, gam[j], beta, k0, -1.0)
        Vp[rows, 2 * j:2 * j + 2] = tangential_fields(am[j], q, gam[j], beta, k0, +1.0)
    vinc = Vp[:, 2 * M:2 * M + 2]                    # order 0, +y, columns s and p

    Zr = np.zeros((4 * N, 2 * N), dtype=complex)
    Asys = np.block([[W * s0[None, :], -Vm, Zr],
                     [(phase_d[:, None] * W) * sd[None, :], Zr, -Vp]])
    rhs = np.vstack([vinc, np.zeros((4 * N, 2), dtype=complex)])
    x = np.linalg.solve(Asys, rhs)
    r = x[4 * N:6 * N].reshape(N, 2, 2)              # [order, out pol, in pol]
    t = x[6 * N:8 * N].reshape(N, 2, 2)

    g0 = gam[M].real
    flux = np.where(gsq > 0.0, gam.real / g0, 0.0)
    DE_R = flux[:, None] * (np.abs(r) ** 2).sum(axis=1)
    DE_T = flux[:, None] * (np.abs(t) ** 2).sum(axis=1)
    k1 = M + 1
    f1 = np.sqrt(gam[k1] / g0) if gsq[k1] > 0.0 else 0.0
    return dict(ms=ms, gam=gam, open=gsq > 0.0, DE_R=DE_R, DE_T=DE_T,
                sum_s=float(DE_R[:, 0].sum() + DE_T[:, 0].sum()),
                sum_p=float(DE_R[:, 1].sum() + DE_T[:, 1].sum()),
                J_T=f1 * t[k1], J_R=f1 * r[k1], t=t, r=r)


def eta_of(J, b):
    """eta(b) = ||J b||^2 for a Jones matrix and a (normalised) incident Jones vector."""
    v = J @ np.asarray(b, dtype=complex)
    return float(np.vdot(v, v).real)


def jones_error(C, J):
    """Report Eq. (63)-(65): the exact worst absolute efficiency error over every normalised
    incident polarisation, e(C, J) = || C^H C - J^H J ||_2."""
    D = C.conj().T @ C - J.conj().T @ J
    D11, D22, D12 = D[0, 0].real, D[1, 1].real, D[0, 1]
    return abs(0.5 * (D11 + D22)) + np.sqrt((0.5 * (D11 - D22)) ** 2 + abs(D12) ** 2)
