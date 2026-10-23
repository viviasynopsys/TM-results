"""rcwa -- rigorous coupled-wave referee, BOTH polarisations, complex index allowed.

Written 2026-10-16 as the polarisation- and loss-aware successor of ``rcwa_te.py``.  The TE
lossless path is the SAME arithmetic as ``rcwa_te.solve`` (same matrix, same anchoring, same
read-out), so every TE table of the 30 September benchmark regenerates unchanged; the TM block
and the complex-index material matrices are those of the separately validated two-polarisation
core (``DBR_volume_gratings/validation/rcw_tm/rcw_rigorous.py``, run 007 of 2026-07-28, and the
explorer package ``psmgrating.rcwa``), transcribed into this benchmark's geometry.
``test_tm_referee.py`` ties the new paths to those two codes and to printed oracle values.

Geometry is F-PSM's (``fpsm.py``): slab 0 <= y <= d, K = K(-sin psi, cos psi) signed, tangential
input wavenumber alpha = beta sin(theta_c - psi).  The field order l carries tangential
wavenumber alpha - l K_x, so F-PSM's signal m = +1 is l = -1.

Material:  n(r) = (n0 - i chi0) + (n1 - i chi1) cos(K.r),  time convention exp(+i w t), so
chi0, chi1 >= 0 absorb.  The exterior is the LOSSLESS real index n0 on both sides (the paper's
Appendix A.1 convention and the benchmark's no-Fresnel-face condition).

    eps_form = "exact"   eps = n^2 with every harmonic (constant n1^2/2, second harmonic n1^2/4)
    eps_form = "linear"  eps = n0^2 + 2 n0 n1 cos(K.r)  (what the OpticStudio RCWA DLL does;
                         diagnostic only, never the reference)

TE (sigma): (h, v) = (E_z, Z0 H_x);  TM (pi): (h, v) = (H_z, E_x / Z0).  TM needs the Toeplitz
matrix of 1/eps as well as of eps; both are built from closed-form Fourier coefficients
(``inv_eps_harmonics``), checked against a 8192-point FFT in the tests.
"""
import numpy as np

__all__ = ["solve", "signal_efficiency", "eps_harmonics", "inv_eps_harmonics", "toeplitz"]

POLS = ("TE", "TM")


def _complex_index(n0, n1, chi0, chi1):
    """Real floats when lossless (keeps the TE path bit-identical to rcwa_te), complex otherwise."""
    n0g = n0 if chi0 == 0.0 else n0 - 1j * chi0
    n1g = n1 if chi1 == 0.0 else n1 - 1j * chi1
    return n0g, n1g


def eps_harmonics(n0g, n1g, eps_form="exact"):
    """{m: eps_m}, m >= 0, of eps(phi) = (n0g + n1g cos phi)^2 (exact) or n0g^2 + 2 n0g n1g cos phi (linear)."""
    if eps_form == "exact":
        return {0: n0g * n0g + n1g * n1g / 2.0, 1: n0g * n1g, 2: n1g * n1g / 4.0}
    if eps_form == "linear":
        return {0: n0g * n0g, 1: n0g * n1g, 2: 0.0}
    raise ValueError(f"eps_form must be 'exact' or 'linear', not {eps_form!r}")


def _branch(a, b):
    """w = sqrt(a^2 - b^2) on the branch for which rho = (w - a)/b has |rho| < 1."""
    w = np.sqrt(complex(a) ** 2 - complex(b) ** 2)
    rho = (w - a) / b
    if abs(rho) >= 1.0:
        w = -w
        rho = (w - a) / b
    return w, rho


def inv_eps_harmonics(n0g, n1g, m, eps_form="exact"):
    """Fourier coefficient g_m of 1/eps(phi), convention 1/eps = sum_m g_m exp(-i m phi).

    exact : 1/(a + b cos)^2 = sum rho^|m| (a + |m| w)/w^3,  w = sqrt(a^2 - b^2), rho = (w - a)/b
    linear: 1/(A + B cos)   = sum rho^|m| / w,              w = sqrt(A^2 - B^2), rho = (w - A)/B
            with A = n0g^2, B = 2 n0g n1g
    Both are even in m, so the exp(-i m phi) and exp(+i m phi) conventions coincide."""
    m = abs(int(m))
    if eps_form == "exact":
        a, b = n0g, n1g
        if b == 0:
            return 1.0 / (a * a) if m == 0 else 0.0
        w, rho = _branch(a, b)
        return rho ** m * (a + m * w) / w ** 3
    if eps_form == "linear":
        A, B = n0g * n0g, 2.0 * n0g * n1g
        if B == 0:
            return 1.0 / A if m == 0 else 0.0
        w, rho = _branch(A, B)
        return rho ** m / w
    raise ValueError(f"eps_form must be 'exact' or 'linear', not {eps_form!r}")


def toeplitz(coeff, N):
    """T[p, l] = c_{p-l} from a callable coeff(m), m = -(N-1)..(N-1)."""
    ms = np.arange(-(N - 1), N)
    cm = np.array([coeff(int(m)) for m in ms], dtype=complex)
    idx = np.arange(N)
    return cm[(idx[:, None] - idx[None, :]) + (N - 1)]


def solve(k0, n0, n1, K, psi, theta_c, d, L=10, eps_form="exact", pol="TE", chi0=0.0, chi1=0.0):
    """One point (scalars; angles in radians, one length unit).  Returns dict with order indices
    ``ls`` and per-order efficiencies ``DE_R``, ``DE_T`` (fractions of the incident power leaving
    through the entry / far face) and ``sum_DE`` (= 1 for a lossless grating)."""
    if pol not in POLS:
        raise ValueError(f"pol must be one of {POLS}, not {pol!r}")
    if chi0 < 0.0 or chi1 < 0.0:
        raise ValueError("chi0, chi1 must be >= 0 (exp(+i w t): a negative imaginary index absorbs)")
    b = k0 * n0
    kx = b * np.sin(theta_c - psi)
    Kx = -K * np.sin(psi)
    Kz = K * np.cos(psi)

    N = 2 * L + 1
    ls = np.arange(-L, L + 1)
    kxl = kx + ls * Kx
    kzl = np.sqrt(b * b - kxl * kxl + 0j)
    kzl = np.where(kzl.imag > 0, -kzl, kzl)                  # Im(k_y,l) <= 0
    kz0 = kzl[L].real

    n0g, n1g = _complex_index(n0, n1, float(chi0), float(chi1))
    eh = eps_harmonics(n0g, n1g, eps_form)
    idx = np.arange(N)
    diff = np.abs(idx[:, None] - idx[None, :])
    E = np.zeros((N, N), dtype=complex)
    for m, v in eh.items():
        E[diff == m] = v

    Lz = np.diag(1j * ls * Kz)
    I = np.eye(N)

    if pol == "TE":
        X2 = np.diag(kxl * kxl)
        M = np.block([[Lz, 1j * k0 * I],
                      [1j * k0 * E - (1j / k0) * X2, Lz]])
        q = kzl / k0
        sgn0, sgnd = -1.0, 1.0
    else:
        X = np.diag(kxl)
        G = toeplitz(lambda m: inv_eps_harmonics(n0g, n1g, m, eps_form), N)
        M = np.block([[Lz, -1j * k0 * E],
                      [(1j / k0) * (X @ G @ X) - 1j * k0 * I, Lz]])
        q = kzl / (k0 * (n0 * n0))                           # real lossless surround
        sgn0, sgnd = 1.0, -1.0

    lam, W = np.linalg.eig(M)
    Wh, Wv = W[:N, :], W[N:, :]

    up = lam.real > 0
    s0 = np.ones(2 * N, dtype=complex)
    sd = np.ones(2 * N, dtype=complex)
    s0[up] = np.exp(-lam[up] * d)
    sd[~up] = np.exp(lam[~up] * d)

    Q = np.diag(q)
    A = np.vstack([(Wv + sgn0 * (Q @ Wh)) * s0[None, :],
                   (Wv + sgnd * (Q @ Wh)) * sd[None, :]])
    rhs = np.zeros(2 * N, dtype=complex)
    rhs[L] = sgn0 * 2.0 * q[L]
    c = np.linalg.solve(A, rhs)

    h0 = Wh @ (c * s0)
    hd = Wh @ (c * sd)
    h0[L] -= 1.0
    DE_R = (kzl.real / kz0) * np.abs(h0) ** 2
    DE_T = (kzl.real / kz0) * np.abs(hd) ** 2
    return dict(ls=ls, DE_R=DE_R, DE_T=DE_T, sum_DE=float(DE_R.sum() + DE_T.sum()))


def signal_efficiency(k0, n0, n1, K, psi, theta_c, d, L=10, eps_form="exact", pol="TE", chi0=0.0, chi1=0.0):
    """Broadcasting wrapper.  Returns (eta_R, eta_T): the l = -1 order's efficiency at the entry
    face and at the far face.  The caller picks the one that matches the geometry."""
    args = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in (k0, n0, n1, K, psi, theta_c, d, chi0, chi1)])
    shape = args[0].shape
    etaR = np.empty(args[0].size)
    etaT = np.empty(args[0].size)
    for i, p in enumerate(zip(*[a.ravel() for a in args])):
        r = solve(*p[:7], L=L, eps_form=eps_form, pol=pol, chi0=float(p[7]), chi1=float(p[8]))
        etaR[i], etaT[i] = r["DE_R"][L - 1], r["DE_T"][L - 1]
    return etaR.reshape(shape), etaT.reshape(shape)
