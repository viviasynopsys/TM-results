# TM reference values for the remote briefs

*Computed 2026-10-16 by `code/reference_values_tm.py` with the gated referee `rcwa.py` (L = 14) and `kogelnik.py`; TE columns reproduce the 30 September benchmark's gate tables. Efficiencies per 1 W incident; angles inside the medium; R+1 / T+1 is the first order on the entry / far face. "eps-linear" is the OpticStudio DLL's permittivity expansion.*


## Converged referee, both polarisations

| case | θ (°) | λ (nm) | face | TE referee | TM referee | TM, eps-linear | TM/TE | TM T0 | Kogelnik TE | Kogelnik TM (κ cos 2θc) |
|---|---|---|---|---|---|---|---|---|---|---|
| `lumerical_R` | 0.0000 | 550 | R+1 | 0.804947 | **0.736870** | 0.736918 | 0.9154 | 0.2629 | 0.805403 | 0.737264 |
| `lumerical_R` | -28.0000 | 550 | R+1 | 0.804947 | **0.736870** | 0.736918 | 0.9154 | 0.2628 | 0.805403 | 0.737264 |
| `lumerical_R` | -14.0000 | 550 | R+1 | 0.023944 | **0.023919** | 0.027389 | 0.9990 | 0.9760 | 0.014014 | 0.014014 |
| `lumerical_R` | 7.5000 | 550 | R+1 | 0.100148 | **0.053714** | 0.053529 | 0.5363 | 0.9460 | 0.084730 | 0.049672 |
| `lumerical_T` | 0.0000 | 550 | T+1 | 0.978518 | **0.916206** | 0.916025 | 0.9363 | 0.0796 | 0.987535 | 0.922215 |
| `lumerical_T` | 5.0000 | 550 | T+1 | 0.013852 | **0.021129** | 0.020901 | 1.5253 | 0.9741 | 0.025587 | 0.029352 |
| `zemax_kba_R` | -45.0000 | 550 | R+1 | 0.873784 | **0.694665** | 0.694597 | 0.7950 | 0.3051 | 0.874565 | 0.695385 |
| `zemax_kba_R` | 0.0000 | 550 | R+1 | 0.873784 | **0.694665** | 0.694597 | 0.7950 | 0.3050 | 0.874565 | 0.695385 |
| `zemax_kba_R` | -22.5000 | 550 | R+1 | 0.012827 | **0.012825** | 0.011817 | 0.9999 | 0.9870 | 0.042365 | 0.042365 |
| `unslanted_T` | 33.8627 | 500 | T+1 | 0.586035 | **0.575192** | 0.575276 | 0.9815 | 0.4247 | 0.585763 | 0.574755 |
| `unslanted_T` | 36.0000 | 500 | T+1 | 0.243709 | **0.025373** | 0.025397 | 0.1041 | 0.9745 | 0.253146 | 0.019329 |
| `unslanted_T_rot5` | 38.8627 | 500 | T+1 | 0.573615 | **0.579853** | 0.579934 | 1.0109 | 0.4200 | 0.573316 | 0.579477 |
| `unslanted_T_rot5` | 41.0000 | 500 | T+1 | 0.243315 | **0.018743** | 0.018821 | 0.0770 | 0.9812 | 0.260150 | 0.014794 |

## The DLL probe of 2026-09-18 (10 layers per period, eps-linear) against the converged TM referee

| θ (°) | DLL TE (10 layers) | DLL TM (10 layers) | DLL TM/TE | referee TM/TE (exact) | referee TM/TE (eps-linear) |
|---|---|---|---|---|---|
| 0.0 | 0.7965196 | 0.7272546 | 0.9130 | 0.9154 | 0.9154 |
| -14.0 | 0.0245959 | 0.0245981 | 1.0001 | 0.9990 | 0.9990 |
| 7.5 | 0.0978717 | 0.0520480 | 0.5318 | 0.5363 | 0.5332 |

## One absorbing point (for the optional Lumerical loss probe)

G1 at Bragg (θ = 0°, 550 nm) with χ₀ = 1×10⁻³ (imaginary part of the mean index, n = 1.5 − 0.001 i), χ₁ = 0:

| pol | first order R+1 | T0 | Σ over all orders (1 − absorbed fraction) |
|---|---|---|---|
| TE | 0.735964 | 0.179353 | 0.915638 |
| TM | 0.668751 | 0.240394 | 0.909354 |

A solver whose imaginary index is silently ignored returns Σ = 1.000000 and the lossless first-order value; that is the gate of the optional step.

