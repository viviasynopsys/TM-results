"""basis -- the geometry and exterior polarisation basis shared by the conical referee and the
conical F-PSM closed forms.  Report: Conical_F_PSM_RCWA_Report.pdf (20 Oct 2026), Sec. 2.A and
Appendix B.A.  Nothing here is physics; it is the convention both codes must share so that their
order-1 Jones matrices are expressed in one basis.

F-PSM coordinates (x', y, z): slab 0 < y < d, exterior index n0 both sides.
    n(xi) = n0 + n1 cos(K xi + phi_g),   xi = -x' sin Psi + y cos Psi
    K_x = K sin Psi,  K_y = K cos Psi,   K = 2 pi / Lambda > 0
    alpha = beta sin theta cos phi,  q = beta sin theta sin phi,  beta = k0 n0        (Eq. 3)
    alpha_m = alpha + m K_x,  gamma_m = sqrt(beta^2 - alpha_m^2 - q^2)              (Eq. 4)
gamma_m > 0 for an open order, -i kappa_m (kappa_m > 0) for a closed one.  Time convention
exp(+i w t), fields exp(-i k.r); sigma = +1 propagates toward +y.

Exterior electric basis of state (m, sigma, rho), report Eq. (92)-(93), bilinear-normalised
(e.e = 1, no conjugation):
    b_m = sqrt(alpha_m^2 + q^2)
    e_s = (q, 0, -alpha_m) / b_m
    e_p = (sigma alpha_m gamma_m, -b_m^2, sigma q gamma_m) / (beta b_m)
    b_m = 0:  e_s = (0, 0, 1),  e_p = e_s x k_hat
"""
import numpy as np

__all__ = ["geometry", "e_vectors", "tangential_fields", "material_harmonics", "B_TOL"]

B_TOL = 1e-12       # relative tolerance below which b_m is treated as the removable singularity


def geometry(k0, n0, K, Psi, theta, phi, ms):
    """Scalars in; per-order arrays out.  Angles in radians, one length unit throughout."""
    beta = k0 * n0
    alpha = beta * np.sin(theta) * np.cos(phi)
    q = beta * np.sin(theta) * np.sin(phi)
    Kx, Ky = K * np.sin(Psi), K * np.cos(Psi)
    ms = np.asarray(ms, dtype=int)
    am = alpha + ms * Kx
    gsq = beta * beta - am * am - q * q
    gam = np.where(gsq >= 0.0, np.sqrt(np.abs(gsq)) + 0j, -1j * np.sqrt(np.abs(gsq)))
    return dict(beta=beta, alpha=alpha, q=q, Kx=Kx, Ky=Ky, ms=ms, am=am, gsq=gsq, gam=gam)


def e_vectors(am, q, gam, beta, sigma):
    """(e_s, e_p) as 3-vectors (x', y, z) for one order; complex when gamma is."""
    b = np.sqrt(am * am + q * q)
    if b <= B_TOL * beta:
        es = np.array([0.0, 0.0, 1.0], dtype=complex)
        khat = np.array([am, sigma * gam, q], dtype=complex) / beta
        ep = np.cross(es, khat)
        return es, ep
    es = np.array([q, 0.0, -am], dtype=complex) / b
    ep = np.array([sigma * am * gam, -b * b, sigma * q * gam], dtype=complex) / (beta * b)
    return es, ep


def tangential_fields(am, q, gam, beta, k0, sigma):
    """For the two polarisations of one exterior plane wave with unit electric amplitude:
    the tangential column (E_x, E_z, g_x, g_z), g = Z0 H = (k x E)/k0, k = (alpha_m, sigma gamma_m, q).
    Returns a 4 x 2 array, columns s and p."""
    es, ep = e_vectors(am, q, gam, beta, sigma)
    k = np.array([am, sigma * gam, q], dtype=complex)
    out = np.empty((4, 2), dtype=complex)
    for j, e in enumerate((es, ep)):
        g = np.cross(k, e) / k0
        out[:, j] = (e[0], e[2], g[0], g[2])
    return out


def material_harmonics(n0, n1, phi_g, hmax):
    """Fourier coefficients of eps = n^2 and of 1/eps for n = n0 + n1 cos(chi + phi_g),
    convention f(chi) = sum_h f_h exp(i h chi), h = -hmax..hmax.  Report Eq. (8)-(9).
    Returns (eps_h, a_h) as dicts keyed by signed h."""
    if abs(n1) >= n0:
        raise ValueError("n1 must be smaller than n0 (the index must stay positive)")
    eps = {0: n0 * n0 + n1 * n1 / 2.0}
    if n1 != 0.0:
        eps[1], eps[-1] = n0 * n1, n0 * n1
        eps[2], eps[-2] = n1 * n1 / 4.0, n1 * n1 / 4.0
    Delta = np.sqrt(n0 * n0 - n1 * n1)
    r = n1 / (n0 + Delta)
    a = {}
    for h in range(-hmax, hmax + 1):
        a[h] = (-r) ** abs(h) * (n0 + abs(h) * Delta) / Delta ** 3
    if phi_g != 0.0:
        eps = {h: v * np.exp(1j * h * phi_g) for h, v in eps.items()}
        a = {h: v * np.exp(1j * h * phi_g) for h, v in a.items()}
    return eps, a
