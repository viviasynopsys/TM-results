# Lumerical RCWA on mixed phase–absorption gratings at conical incidence — verdict

*2026-10-23, solver PC, Ansys Lumerical 2026 R1, build 8.35.4609 (lumapi from `C:\Program Files\ANSYS Inc\v261\Lumerical\api\python`), Python 3.13.14 (numpy 2.5.1, scipy 1.18.1, matplotlib 3.11.1), Windows. British English.*

## Verdict

**Lumerical's RCWA solver handles a complex (absorbing) index through an `import` primitive, and with it agrees with the independent complex-permittivity conical referee on all 58 configurations of sets A, B and C: efficiencies and absorbed fractions to ≤ 2.1 × 10⁻⁴ rms (≤ 5.1 × 10⁻⁵ on the 54 grating configurations of sets A and B), every power sum below 1, every Jones comparison within the triangle bound, and the author's 30 certified bounds E19 confirmed. VERIFIED 2026-10-23.**

**The spatial index equation does not accept a complex index: REJECTED, VERIFIED 2026-10-23.** Gates X0–X4 and X6 pass. X5 (set D phase scan) was not run: UNTESTED.

## Which route, sign and token worked

| item | result |
|---|---|
| route `equation` | **REJECTED.** For every token (`1i`, `i`, `1j`) and both signs the property reads back as the real equation `1.5+0*cos(...)`; the imaginary part is never stored (`results/lumerical_mixed_probe_loss_equation.FAILED.txt`) |
| route `import` | **WORKS.** `addimport` + `importnk2(n, x, y, z)` on a 256 × 2 × 1024 grid, mesh order properties set without error, rectangle index set to the plain n0 |
| sign | Lumerical's index is `n + i κ` (k sign +1): absorbs, power sums below 1. VERIFIED 2026-10-23 |
| token | `n+ik array` (the import route has no imaginary literal) |
| P1 | uniform slab n0 − i 0.02, 6 µm, `lumerical_R` model: max \|T0 − exact Fresnel oracle\| = **1.44 × 10⁻¹³** (gate 10⁻³). T0 θ = 0: 0.064463 (s), 0.064463 (p); θ = 20°, φ = 150°: 0.054065 / 0.054063; power sums 0.064508 / 0.064508 and 0.054122 / 0.054096, all equal to the brief's table to six figures |
| P2 | published absorbing point, θ = 0, 160-layer staircase: R+1 TE **0.735904**, TM **0.668684**, power sums **0.915632 / 0.909348** against the referee 0.735964 / 0.668751, 0.915638 / 0.909354. Max deviation **6.72 × 10⁻⁵** (gate 2 × 10⁻³). PASS |

`results/loss_route.json` records route `import`, k sign +1, harmonics 3, layers 160.

## Gates, from `results/summary_mixed.csv` (58 configs: A 30, B 24, C 4)

| gate | requirement | result |
|---|---|---|
| X0 | `power_max` < 1; absorbed-fraction rms ≤ 10⁻³ | **PASS.** `power_max` ≤ 0.9293 (C03); absorbed rms ≤ 4.6 × 10⁻⁵ (A16), max absolute absorbed error 1.1 × 10⁻⁴ (set C) |
| X1 | `rms_s`, `rms_p` ≤ 10⁻³ | **PASS.** Set A ≤ 4.3 × 10⁻⁵ (s), 5.1 × 10⁻⁵ (p); set B ≤ 7 × 10⁻⁶, 9 × 10⁻⁶; set C ≤ 1.84 × 10⁻⁴, 2.10 × 10⁻⁴ (worst: C02, the near-grazing solver-vs-referee test) |
| X2 | amplitude convention resolved; time convention; Jones error ≤ 10⁻³ | **PASS.** Amplitude convention "in-first, power-normalised" on all 58. Time convention `exp(−iωt): amplitudes conjugated` on 57; the single `exp(+iωt)` entry (A30) and the two `p_source_sign` −1 entries (A10, A20) are the zero-signal one-sided configs, where the signal is 10⁻¹³–10⁻¹⁶ and no convention can be decided. `jones_err_lumerical_vs_referee` max **3.55 × 10⁻⁴** (C02) |
| X3 | `jones_within_bound` | **PASS**, True on 58 of 58 |
| X4 | set A `within_E19`; Lumerical column within a factor two | **PASS**, True on 30 of 30; corrected error at θ0 / E19 in [0.99984, 1.00000141] over the non-zero bounds; (corrected-vs-Lumerical) / (corrected-vs-referee + Lumerical-vs-referee) in [0.885, 1.000] |
| X5 | set D phase scan | **UNTESTED** (set D not run) |
| X6 | convergence | **PASS**, see below |

Per-set worst values:

| set | n | rms_s | rms_p | absorbed rms s / p | Jones err | `power_max` |
|---|---|---|---|---|---|---|
| A | 30 | 4.3 × 10⁻⁵ | 5.1 × 10⁻⁵ | 3.9 × 10⁻⁵ / 4.6 × 10⁻⁵ | 5.5 × 10⁻⁵ | 0.523 |
| B | 24 | 7 × 10⁻⁶ | 9 × 10⁻⁶ | 8 × 10⁻⁶ / 1.0 × 10⁻⁵ | 1.1 × 10⁻⁵ | 0.050 |
| C | 4 | 1.84 × 10⁻⁴ | 2.10 × 10⁻⁴ | 2.2 × 10⁻⁵ / 2.4 × 10⁻⁵ | 3.55 × 10⁻⁴ | 0.929 |

One-sided configs (δ = 270°, κ1 = n1): signal efficiency from Lumerical A10 ≤ 9.8 × 10⁻¹⁵, A20 ≤ 1.8 × 10⁻¹⁶, A30 ≤ 5.5 × 10⁻¹³ (brief: ≤ 10⁻⁶); the referee gives zero. Their absorbed fractions agree with the referee to 1.1 × 10⁻⁸, 1.3 × 10⁻⁷ and 1.5 × 10⁻⁷.

Per-config numbers are in `results/summary_mixed.csv`; all 58 figures `figures/mixed_<case_id>.png` were looked at: Lumerical markers lie on the referee curves in every panel; the corrected closed form departs only at the near-grazing edge of C02, as the brief says.

## Convergence ladder (X6), A05, θ0 = 25°, φ = 150°

| harmonics, layers | η_s | η_p | absorbed s | absorbed p |
|---|---|---|---|---|
| referee (21 orders) | 0.0489532 | 0.0336021 | 0.9053685 | 0.9124841 |
| 3, 80 | 0.0489301 | 0.0335852 | 0.9053881 | 0.9124966 |
| 3, 160 | 0.0489448 | 0.0335959 | 0.9053757 | 0.9124888 |
| 3, 320 | 0.0489483 | 0.0335985 | 0.9053724 | 0.9124866 |
| 5, 160 | 0.0489448 | 0.0335959 | 0.9053757 | 0.9124888 |

Harmonics 3 → 5 change nothing (to the printed 7 figures). Layers 80 → 160 → 320: η_s steps 1.47 × 10⁻⁵ then 3.5 × 10⁻⁶ (ratio 4.2), monotonic, the 1/N² staircase of the lossless leg. 3 harmonics and 160 layers kept.

## δ orientation (`results/delta_orientation.json`, VERIFIED from the data 2026-10-23)

| geometry (φ = 150°) | Lumerical one-sided η_s | referee with +δ / with −δ | decision |
|---|---|---|---|
| slanted T | 0.037511 | 4.8 × 10⁻²⁸ / 0.037520 | δ sign −1 |
| slanted R | 9.0 × 10⁻¹⁷ | 2.6 × 10⁻²⁷ / 0.040528 | δ sign +1 |
| unslanted T | 0.038745 | 1.1 × 10⁻²⁷ / 0.038749 | δ sign −1 |

θ sign +1 and signal order +1 on all three, as on every earlier grating. The fourth geometry named in the brief (unslanted T at φ = 135°) belongs to set D and was not needed.

## What was changed, and why

1. `code/lumerical_rcwa_mixed.py`, `remove_import()`: added `fdtd.switchtolayout()` before the import object is selected and deleted. The first import-route probe stopped with `you cannot modify simulation objects while in analysis mode, use switchtolayout first` after P01 had passed, because the session is in analysis mode after a solve. Nothing else in the logic changed; P01 returned identical numbers before and after.
2. The same file's docstring and two inline labels now state the results above with the date; the `UNTESTED`/`ASSUMED` labels on the complex index, the equation tokens and `importnk2` are replaced by the evidence. Nothing else (`lumerical_rcwa_conical.py`, `lumerical_rcwa.py`) was touched.
3. No gate was loosened, no number tuned, no `.fsp` hand-edited, no file under `reference/`, `lumerical/ansys_example/`, `code/david_mixed/`, `evidence/` or the two `reference_*.csv` changed.

Step 0 on this PC: `check_env.py lumerical` READY; `test_conical.py` 0 failures of 12 gates (1287 s); `test_mixed.py` 0 failures of 11 gates (665 s); `--list` 100 configs; `compare_mixed.py` "no Lumerical table yet".

## Observations

- **Unslanted gratings are not at the brief's 10⁻¹¹.** The brief expected about 10⁻¹¹ on the unslanted configs. Measured: 10⁻⁷ to 10⁻⁵ (A21–A29 up to 6.4 × 10⁻⁶ rms; B17–B24 ≤ 6.7 × 10⁻⁷; C04, the K_z = 0 bypass, 1.4 × 10⁻⁵). All pass the 10⁻³ gate by two orders of magnitude. The cause is UNTESTED; the likely one is the 256-point sampling of the profile in the import grid (the lossless unslanted runs used the exact equation), but the grid was not varied.
- **Set C02** has its Bragg cone at 69.7°; it tests the solver against the referee only, where it is the loosest of the 58 (2.1 × 10⁻⁴ rms) and still passes.
- **The brief said the equation route might be ignored; it is worse than ignored:** the equation text is not stored at all, so the probe's read-back check fails on the first call rather than the power sum coming out 1.

## What remains UNTESTED

- Set D (40 configs): D01 θ sweep and the X5 relative-phase scan, hence the δ-orientation of unslanted T at φ = 135°.
- Whether the import grid (256 × 2 × 1024) is converged: only the layer and harmonic ladders were run (X6).
- Lumerical licence edition: not read out; the RCWA solver ran with no licence dialog or error on this PC. Edition string UNTESTED.
- Any index with κ larger than 0.02 (P01) or a thickness above 13.3 µm.
- `p_source_sign` and the time convention on non-zero signal beyond what the 55 non-null configs give: ASSUMED +1 and `exp(−iωt)` as found.
