"""kogelnik_te -- the 30 September 2026 Kogelnik engine, lossless TE, kept verbatim as the regression oracle for kogelnik.py (2026-10-16).

H. Kogelnik, Bell Syst. Tech. J. 48, 2909 (1969), eqs. (42)-(43) [transmission]
and (57)-(58) [reflection], written in F-PSM's geometry (``fpsm.py``) so the
three models take identical arguments:

    rho = beta (sin th, cos th),   th = theta_c - psi   (input wave, from the normal y)
    K   = K (-sin psi, cos psi)    signed
    sigma = rho - K                (Kogelnik's forced signal wave)
    cR = cos th,  cS = cos th - (K / beta) cos psi
    dephasing  vartheta = (beta^2 - |sigma|^2) / (2 beta)
    kappa = pi n1 / lam = k0 n1 / 2

This is the method OpticStudio's native volume-hologram model and the Lumerical
article's "Kogelnik" curve use.  cS < 0 -> reflection, cS > 0 -> transmission.
"""
import numpy as np

__all__ = ["efficiency"]


def efficiency(k0, n0, n1, K, psi, theta_c, d):
    k0, n0, n1, K, psi, theta_c, d = np.broadcast_arrays(
        *[np.asarray(v, dtype=float) for v in (k0, n0, n1, K, psi, theta_c, d)])
    beta = k0 * n0
    th = theta_c - psi
    sx = beta * np.sin(th) + K * np.sin(psi)
    sy = beta * np.cos(th) - K * np.cos(psi)
    cR = np.cos(th)
    cS = sy / beta
    vth = (beta ** 2 - (sx ** 2 + sy ** 2)) / (2.0 * beta)
    kappa = 0.5 * k0 * n1

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        xi = vth * d / (2.0 * cS)
        # transmission, eqs. (42)-(43)
        nuT = kappa * d / np.sqrt(np.abs(cR * cS))
        rootT = np.sqrt(nuT ** 2 + xi ** 2)
        etaT = np.where(rootT == 0, 0.0, np.sin(rootT) ** 2 / (1.0 + xi ** 2 / nuT ** 2))
        # reflection, eqs. (57)-(58): nu is imaginary, written with |cS|
        nuR = nuT
        rootR = np.sqrt(nuR ** 2 - xi ** 2 + 0j)
        sh2 = np.sinh(rootR) ** 2
        etaR = np.real(1.0 / (1.0 + (1.0 - xi ** 2 / nuR ** 2) / sh2))
    eta = np.where(cS < 0, etaR, etaT)
    # Kogelnik's forced wave sigma knows nothing of cut-off: past it the formula
    # still returns a number although no first order can propagate.  OpticStudio's
    # native model returns 0 there (VERIFIED 2026-09-17, zemax_kba_R theta >= 17.5
    # deg); do the same so the curves mean the same thing.
    open_signal = beta ** 2 - sx ** 2 > 0.0
    return np.where((n1 == 0) | ~open_signal, 0.0, eta)
