"""rcwa_te -- rigorous coupled-wave reference, TE (sigma), lossless, index-matched.

The referee for F-PSM.  TE-only reduction of a separately validated two-polarisation RCWA core: slab
0 <= y <= d with n = n0 + n1 cos(K.r), exterior index n0 on both sides, Floquet
orders |l| <= L, solved by diagonalisation with eigenmode anchoring (no
growing exponentials).  ``test_engines.py`` checks it against the explorer
package (when importable) and against energy conservation.

Geometry is F-PSM's (``fpsm.py``): K = K(-sin psi, cos psi) signed, tangential
input wavenumber alpha = beta sin(theta_c - psi).  The field order l here
carries tangential wavenumber alpha - l K_x, so F-PSM's signal m = +1 is
l = -1.
"""
import numpy as np

__all__ = ["solve", "signal_efficiency"]


def solve(k0, n0, n1, K, psi, theta_c, d, L=10, eps_form="exact"):
    """One point (scalars; angles in radians, one length unit).  Returns dict
    with order indices ``ls`` and per-order efficiencies ``DE_R``, ``DE_T``.

    ``eps_form``: "exact" (default, unchanged behaviour) expands eps = (n0 + n1 cos)^2 with its
    DC term n1^2/2 and second harmonic n1^2/4.  "linear" uses eps = n0^2 + 2 n0 n1 cos, i.e. drops
    every n1^2 term -- this is what the Ansys Zemax RCWA hologram DLL does (VERIFIED 2026-09-18:
    with it, this code reproduces vhg_constant_hologram230313.dll, converged in z-slices, to
    ~1e-5 at Bragg and at the sidelobe minimum where the two forms differ by 3.5e-3).  It exists
    only to diagnose that difference; the referee is the exact form."""
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

    # eps = (n0 + n1 cos)^2 : harmonics 0, +-1, +-2
    if eps_form == "exact":
        eh = {0: n0 * n0 + n1 * n1 / 2.0, 1: n0 * n1, 2: n1 * n1 / 4.0}
    elif eps_form == "linear":
        eh = {0: n0 * n0, 1: n0 * n1, 2: 0.0}
    else:
        raise ValueError(f"eps_form must be 'exact' or 'linear', not {eps_form!r}")
    idx = np.arange(N)
    diff = np.abs(idx[:, None] - idx[None, :])
    E = np.zeros((N, N), dtype=complex)
    for m, v in eh.items():
        E[diff == m] = v

    Lz = np.diag(1j * ls * Kz)
    X2 = np.diag(kxl * kxl)
    I = np.eye(N)
    M = np.block([[Lz, 1j * k0 * I],
                  [1j * k0 * E - (1j / k0) * X2, Lz]])
    lam, W = np.linalg.eig(M)
    Wh, Wv = W[:N, :], W[N:, :]

    up = lam.real > 0
    s0 = np.ones(2 * N, dtype=complex)
    sd = np.ones(2 * N, dtype=complex)
    s0[up] = np.exp(-lam[up] * d)
    sd[~up] = np.exp(lam[~up] * d)

    Q = np.diag(kzl / k0)
    A = np.vstack([(Wv - Q @ Wh) * s0[None, :],
                   (Wv + Q @ Wh) * sd[None, :]])
    rhs = np.zeros(2 * N, dtype=complex)
    rhs[L] = -2.0 * kzl[L] / k0
    c = np.linalg.solve(A, rhs)

    h0 = Wh @ (c * s0)
    hd = Wh @ (c * sd)
    h0[L] -= 1.0
    DE_R = (kzl.real / kz0) * np.abs(h0) ** 2
    DE_T = (kzl.real / kz0) * np.abs(hd) ** 2
    return dict(ls=ls, DE_R=DE_R, DE_T=DE_T, sum_DE=float(DE_R.sum() + DE_T.sum()))


def signal_efficiency(k0, n0, n1, K, psi, theta_c, d, L=10, eps_form="exact"):
    """Broadcasting wrapper.  Returns (eta_R, eta_T): the l = -1 order's
    efficiency at the entry face and at the far face.  The caller picks the
    one that matches the geometry (F-PSM's case flag)."""
    args = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in (k0, n0, n1, K, psi, theta_c, d)])
    shape = args[0].shape
    etaR = np.empty(args[0].size)
    etaT = np.empty(args[0].size)
    for i, p in enumerate(zip(*[a.ravel() for a in args])):
        r = solve(*p, L=L, eps_form=eps_form)
        etaR[i], etaT[i] = r["DE_R"][L - 1], r["DE_T"][L - 1]
    return etaR.reshape(shape), etaT.reshape(shape)
