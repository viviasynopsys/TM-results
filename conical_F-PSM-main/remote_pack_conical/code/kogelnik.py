"""kogelnik -- Kogelnik's 1969 two-wave coupled-wave efficiency, TE or TM, lossless or absorbing.

H. Kogelnik, Bell Syst. Tech. J. 48, 2909 (1969).  Written in F-PSM's geometry (``fpsm.py``) so
the models take identical arguments:

    rho = beta (sin th, cos th),   th = theta_c - psi   (input wave, from the normal y)
    K   = K (-sin psi, cos psi)    signed
    sigma = rho - K                (Kogelnik's forced signal wave)
    cR = cos th,  cS = cos th - (K / beta) cos psi
    dephasing  vartheta = (beta^2 - |sigma|^2) / (2 beta)
    kappa = (k0 / 2) (n1 - i chi1)              complex for an absorption-modulated grating
    alpha = k0 chi0                             amplitude absorption constant of the mean index
    TM (pi, "H-mode"): kappa -> kappa cos(2 theta_c), theta_c the angle between the input wave and K

The coupled equations (Kogelnik eqs. 19-20 with loss, time convention exp(+i w t)):

    cR R' + alpha R            = -i kappa S
    cS S' + (alpha + i vth) S  = -i kappa R

Boundary conditions R(0) = 1 and S(0) = 0 (transmission, cS > 0) or S(d) = 0 (reflection,
cS < 0); eta = |cS / cR| |S|^2 at the exit face.  Two code paths:

  * ``pol="TE"`` with chi0 = chi1 = 0 runs the 30 September 2026 lossless formulae UNCHANGED
    (Kogelnik eqs. 42-43 and 57-58), so every TE table of the benchmark regenerates bit for bit
    (``kogelnik_te.py`` is that engine kept verbatim; ``test_kogelnik_pol.py`` checks equality).
  * everything else solves the 2 x 2 constant-coefficient system in closed form (its matrix
    exponential), which reduces to the same formulae when the loss vanishes -- checked to 1e-12
    against the TE path and to 1e-10 against the explorer package's Kogelnik (TM, lossy) by
    ``test_kogelnik_pol.py``.

This is the model OpticStudio's native volume-hologram object and the Lumerical article's
"Kogelnik" curve use (TE; whether OpticStudio applies the TM factor is UNTESTED 2026-10-16).
Past the first order's cut-off (sigma evanescent) the efficiency is reported as 0, as
OpticStudio's native model does (VERIFIED 2026-09-17).
"""
import numpy as np

__all__ = ["efficiency"]

POLS = ("TE", "TM")


def _geometry(k0, n0, K, psi, theta_c):
    beta = k0 * n0
    th = theta_c - psi
    sx = beta * np.sin(th) + K * np.sin(psi)
    sy = beta * np.cos(th) - K * np.cos(psi)
    cR = np.cos(th)
    cS = sy / beta
    vth = (beta ** 2 - (sx ** 2 + sy ** 2)) / (2.0 * beta)
    return beta, sx, cR, cS, vth


def _lossless_te(k0, n0, n1, K, psi, theta_c, d):
    """The 30 September engine, verbatim (kogelnik_te.efficiency)."""
    beta, sx, cR, cS, vth = _geometry(k0, n0, K, psi, theta_c)
    kappa = 0.5 * k0 * n1
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        xi = vth * d / (2.0 * cS)
        nuT = kappa * d / np.sqrt(np.abs(cR * cS))
        rootT = np.sqrt(nuT ** 2 + xi ** 2)
        etaT = np.where(rootT == 0, 0.0, np.sin(rootT) ** 2 / (1.0 + xi ** 2 / nuT ** 2))
        nuR = nuT
        rootR = np.sqrt(nuR ** 2 - xi ** 2 + 0j)
        sh2 = np.sinh(rootR) ** 2
        etaR = np.real(1.0 / (1.0 + (1.0 - xi ** 2 / nuR ** 2) / sh2))
    eta = np.where(cS < 0, etaR, etaT)
    open_signal = beta ** 2 - sx ** 2 > 0.0
    return np.where((n1 == 0) | ~open_signal, 0.0, eta)


def _sinhc(x):
    """sinh(x)/x, with its series below |x| = 1e-6 (complex x)."""
    small = np.abs(x) < 1e-6
    xs = np.where(small, 1.0, x)
    return np.where(small, 1.0 + x * x / 6.0, np.sinh(xs) / xs)


def _general(k0, n0, n1, K, psi, theta_c, d, pol, chi0, chi1):
    """Closed-form solution of the 2 x 2 system for any polarisation and loss.

    With M = [[-alpha/cR, -i kappa/cR], [-i kappa/cS, -(alpha + i vth)/cS]] and tau = tr(M)/2,
    Delta = sqrt((m11 - m22)^2/4 + m12 m21):  exp(M d) = e^{tau d} [cosh(Delta d) I + sinh(Delta d)/Delta (M - tau I)].
    Transmission: S(d) = E21.  Reflection: S(d) = E21 + E22 S(0) = 0 -> S(0) = -E21/E22."""
    beta, sx, cR, cS, vth = _geometry(k0, n0, K, psi, theta_c)
    kappa = 0.5 * k0 * (n1 - 1j * chi1)
    if pol == "TM":
        kappa = kappa * np.cos(2.0 * theta_c)
    alpha = k0 * chi0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        m11 = -alpha / cR + 0j
        m12 = -1j * kappa / cR
        m21 = -1j * kappa / cS
        m22 = -(alpha + 1j * vth) / cS
        tau = 0.5 * (m11 + m22)
        Delta = np.sqrt(0.25 * (m11 - m22) ** 2 + m12 * m21)
        x = Delta * d
        # exp(tau d) can be large for a reflection grating (cS < 0 makes -alpha/cS > 0); the ratio
        # E21/E22 is what matters there, so factor e^{tau d} out and keep only the bracket B.
        sh = _sinhc(x) * d                       # sinh(Delta d)/Delta
        ch = np.cosh(x)
        B21 = sh * m21                           # bracket element (2,1): (M - tau I)_21 = m21
        B22 = ch + sh * (m22 - tau)
        etaT = np.abs(cS / cR) * np.abs(np.exp(tau * d) * B21) ** 2
        etaR = np.abs(cS / cR) * np.abs(B21 / B22) ** 2
    eta = np.where(cS < 0, etaR, etaT)
    open_signal = beta ** 2 - sx ** 2 > 0.0
    return np.where(((n1 == 0) & (chi1 == 0)) | ~open_signal, 0.0, np.real(eta))


def efficiency(k0, n0, n1, K, psi, theta_c, d, pol="TE", chi0=0.0, chi1=0.0):
    """First-order efficiency, 0 past the signal's cut-off.  Angles in radians, one length unit;
    chi0, chi1 the imaginary parts of the mean index and of the modulation (>= 0 absorb)."""
    if pol not in POLS:
        raise ValueError(f"pol must be one of {POLS}, not {pol!r}")
    k0, n0, n1, K, psi, theta_c, d, chi0, chi1 = np.broadcast_arrays(
        *[np.asarray(v, dtype=float) for v in (k0, n0, n1, K, psi, theta_c, d, chi0, chi1)])
    if np.any(chi0 < 0) or np.any(chi1 < 0):
        raise ValueError("chi0, chi1 must be >= 0")
    if pol == "TE" and not np.any(chi0) and not np.any(chi1):
        return _lossless_te(k0, n0, n1, K, psi, theta_c, d)
    return _general(k0, n0, n1, K, psi, theta_c, d, pol, chi0, chi1)
