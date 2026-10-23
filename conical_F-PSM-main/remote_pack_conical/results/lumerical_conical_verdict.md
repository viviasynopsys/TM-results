# Lumerical RCWA at conical incidence: verdict

*Written 2026-10-23 on the work PC. Lumerical 8.35.4609 (Ansys v261), Python 3.13.14, numpy 2.5.1, scipy 1.18.1, matplotlib 3.11.1. British English.*

**VERDICT. Lumerical's RCWA solver handles conical incidence correctly. On all seven gratings, 14 (case, φ) sweeps and 601 solved rows, its efficiencies agree with the independent conical referee to rms ≤ 4.9 × 10⁻⁵ (≤ 2.1 × 10⁻⁴ on the grazing `lumerical_T` row), power is conserved to 10⁻⁶ in every sweep (10⁻¹² at the Bragg points), and the s→p conversion of the conical signature is reproduced (VERIFIED 2026-10-23). Gates C0, C1, C2, C3 and C4 pass. Gate C5 is only partly met: the amplitude convention resolves (in-first, power-normalised) on all 14 rows, but the Jones error against the referee exceeds 10⁻³ on 7 of 14 rows, and the corrected-model Jones comparison passes on only 3 of 14 (the two planar φ = 180° rows and `lumerical_T`). The cause is a relative phase in the s–p coherence of Lumerical's amplitudes that a ±1 source sign does not remove; no normalisation was fitted or guessed. Details and what to establish next are below.**

## Gate results (VERIFIED 2026-10-23, from `results/summary_conical.csv`)

| gate | what | result |
|---|---|---|
| Step 0 | `check_env.py lumerical` / `test_conical.py` / `test_fpsm_vectors.py` / `compare_conical.py lumerical_R` | READY / 0 failures of 12 gates (362 s) / 0 failures / "no Lumerical table yet" |
| C1 | planar fit | sign +1 on every case; signal order +1 (−1 on BENCH `unslanted_T`); rms 1.6 × 10⁻⁵ to 4.4 × 10⁻⁵ (1.3 × 10⁻¹² and 8.7 × 10⁻¹² on the two unslanted gratings) |
| C0 | power sums, s and p | 1.000000 to six figures on every row of every sweep; ≤ 3 × 10⁻¹³ from 1 at the Bragg points |
| C3 | conical signature, `lumerical_R`, φ = 150°, θ = 24.369° | η(s in) = 0.787212, of which **s→p 0.787177** (s→s 3.6 × 10⁻⁵); η(p in) = 0.753131, **p→s 0.751594** (referee 0.787246 / 0.787211 and 0.753167 / 0.751631). A solver ignoring φ would give 0.805 with no conversion. φ = 100°: s→p 0.738752, p→s 0.802420 |
| C4 | convergence (below) | harmonics 3 → 5 changes ≤ 4.1 × 10⁻¹⁰; layers 80 → 160 → 320 monotone from below with difference ratio 4.00 on all 8 values |
| C2 | rms s, p vs referee | 14 of 14 pass (≤ 1 × 10⁻³) |
| C5 | amplitude convention resolved | 14 of 14 ("in-first, power-normalised") |
| C5 | Jones error Lumerical vs referee ≤ 10⁻³ | 7 of 14 pass; **7 fail** |
| C5 | corrected-model Jones error within a factor two of `corr_max_err` | 3 of 14 pass; 11 fail |

| case | φ | n | rms s | rms p | corrected model err | Jones L–R | Jones C–L | C2 | C0 | C5 (Jones) |
|---|---|---|---|---|---|---|---|---|---|---|
| lumerical_R | 100 | 42 | 3.82e-5 | 3.43e-5 | 8.5e-5 | 2.48e-4 | 2.46e-4 | pass | pass | L–R pass, C–L fail |
| lumerical_R | 120 | 49 | 4.04e-5 | 3.86e-5 | 1.9e-5 | 5.76e-4 | 5.83e-4 | pass | pass | L–R pass, C–L fail |
| lumerical_R | 150 | 49 | 4.40e-5 | 4.28e-5 | 2.8e-5 | 5.18e-4 | 5.38e-4 | pass | pass | L–R pass, C–L fail |
| lumerical_R | 180 | 49 | 4.52e-5 | 4.34e-5 | 1.4e-4 | 5.99e-5 | 1.18e-4 | pass | pass | pass |
| lumerical_T | 100 | 49 | 1.89e-4 | 2.09e-4 | 5.75e-1 (outside the closed forms' domain, expected) | 5.70e-2 | 5.76e-1 | pass | pass | L–R fail, C–L pass |
| zemax_kba_R | 100 | 49 | 4.00e-5 | 2.81e-5 | 1.6e-4 | 1.72e-3 | 1.75e-3 | pass | pass | fail |
| zemax_kba_R | 120 | 49 | 4.18e-5 | 3.82e-5 | 6.0e-4 | 3.94e-3 | 3.95e-3 | pass | pass | fail |
| zemax_kba_R | 150 | 49 | 4.55e-5 | 4.19e-5 | 3.7e-5 | 3.43e-3 | 3.45e-3 | pass | pass | fail |
| zemax_kba_R | 180 | 49 | 4.82e-5 | 4.00e-5 | 1.7e-4 | 7.56e-5 | 1.91e-4 | pass | pass | pass |
| unslanted_T | 30 | 49 | 5.4e-12 | 9.9e-12 | 4.7e-6 | 1.34e-2 | 1.34e-2 | pass | pass | fail |
| unslanted_T | 45 | 49 | 2.1e-11 | 8.5e-12 | 5.7e-6 | 2.03e-2 | 2.03e-2 | pass | pass | fail |
| Slanted_T | 150 | 23 | 2.70e-5 | 2.27e-5 | 1.4e-6 | 1.97e-3 | 1.97e-3 | pass | pass | fail |
| Slanted_R | 150 | 23 | 4.14e-5 | 4.14e-5 | 1.5e-6 | 2.59e-4 | 2.58e-4 | pass | pass | L–R pass, C–L fail |
| Unslanted_T (DAVID) | 150 | 23 | 8.3e-12 | 5.9e-12 | 9.2e-7 | 5.59e-4 | 5.59e-4 | pass | pass | L–R pass, C–L fail |

Bragg-point efficiencies against the brief's table: every row within 4.7 × 10⁻⁵ of the referee except `lumerical_T` (1.8 × 10⁻⁴ in s, 1.0 × 10⁻⁴ in p; near-grazing signal order, γ₁/β = 0.22). The slanted gratings sit about 3.5 × 10⁻⁵ below the referee, which is the 160-layer staircase: the ladder moves them up by 2.8 × 10⁻⁵ at 320 layers. Off-Bragg discriminators: `lumerical_R`, φ = 150°, θ = 26.369°: 0.758266 / 0.700240 (referee 0.758304 / 0.700280); `unslanted_T`, φ = 30°, θ = 42.046°: 0.141704 / 0.084094 (referee 0.141704 / 0.084094). The φ = 180° rows reproduce the September planar values (`lumerical_R` 0.804914 / 0.736832 against 0.804947 / 0.736870; `zemax_kba_R` 0.873758 / 0.694626 against 0.873784 / 0.694665), so the continuity check with the September pack holds (VERIFIED 2026-10-23).

## Convergence ladder (`results/lumerical_conical_converge.txt`)

`lumerical_R`, Bragg point of each azimuth. Layers 80 → 160 → 320 at 3 harmonics: the change is +9.8 × 10⁻⁵ to +1.1 × 10⁻⁴ and then +2.5 × 10⁻⁵ to +2.8 × 10⁻⁵ (ratio 4.00 on all eight s and p values); harmonics 3 → 5 at 160 layers changes nothing above 4.1 × 10⁻¹⁰. Settings kept: 3 harmonics, 160 layers (VERIFIED 2026-10-23 for `lumerical_R`; the other six cases are checked at 160 layers only through their agreement with the referee).

## Amplitude convention and the C5 shortfall

Found (VERIFIED 2026-10-23, by `compare_conical.py`): `grating_characterization` amplitudes are **power-normalised** (|amplitude|² is the order's power) and **in-first**: `Rsp`/`Tsp` is incident s, outgoing p. The out-first reading was tested as well on the three gratings where s→p equals p→s (`unslanted_T`, DAVID `Unslanted_T`) and where the power test cannot tell the readings apart: it breaks the diagonal terms of J†J against the referee (for example 0.016567 against 0.014513), so in-first stands.

The shortfall: the Jones error compares J†J with the referee's. Its diagonals (η_s, η_p) and the **magnitude** of the s–p coherence |J†J₀₁| match to ≤ 2.5 × 10⁻⁵ (2 × 10⁻¹¹ on the unslanted gratings), but its **phase** does not. Measured relative phase of the coherence, Lumerical against referee: `lumerical_R` φ = 150°: −0.17° to −0.43°; `zemax_kba_R` φ = 120°: −1.0° to −1.2°; `Slanted_T`: about −0.4°; `unslanted_T` φ = 30°: −0.5° to −2.2° and **+6.3° at Bragg**. It varies with θ, so it is not the ±1 `p_source_sign` that `compare_conical.py` tries (reported as +1 on every row; ASSUMED). The planar φ = 180° rows, which have no s–p coherence, pass. I have not fitted or guessed a normalisation. Hypotheses, none tested: a θ-dependent phase in Lumerical's p-source definition, or a phase reference of the incident or diffracted order that differs between Lumerical and the referee's basis. To establish on the benchmark PC: the first ten rows of `grating_characterization` for `unslanted_T` φ = 30° are in `results/lumerical_conical_unslanted_T.csv` (columns `Rss_re` … `Tpp_im`).

## What was changed in the driver and why (all in `code/`)

The September driver `lumerical_rcwa.py` is unchanged, and no entry of the NAME MAP `P` or `AMP_KEYS` was corrected: every name the brief lists was present on this build (probe, step 1).

1. **`lumerical_rcwa_conical.py`, `planar_sign_fit`**: for the DAVID gratings the planar fit sampled θ = θ_conical ± 3° at φ_lum = 0, where `Slanted_T` has no Bragg point (its planar Bragg is at Lumerical θ = −21.19°, equal to the recording beam). The first `Slanted_T` run therefore fitted order −1 at rms 2.9 × 10⁻³ and failed with "signal order −1 not in this solve's order axis [0, 1]" (first-attempt table was empty, `*.FAILED.csv`). Diagnosed by a planar scan: at θ = −21.19°, order +1 gives 0.81612 against the referee's 0.81618. The fit now centres on the planar Bragg angle, signed in Lumerical's axes. After the change the planar fit gives sign +1, order +1, rms 2.6 × 10⁻⁵ (`Slanted_T`), 4.4 × 10⁻⁵ (`Slanted_R`).
2. **`lumerical_rcwa_conical.py`, `configure_conical`**: the K_z = 0 bypass was selected by the name `unslanted_T` only, so the DAVID `Unslanted_T` (Ψ = 90°) went down the setup-script path and failed with "asked 160 layers, solver made 41 interfaces". The bypass is now also used whenever cos Ψ ≈ 0, with the DAVID recording beams. It is the same verified September routine; the recording angle is read from the family's own `MODEL` entry in memory only.
3. **`cases_conical.py` (`file_stem`), `lumerical_rcwa_conical.py`, `compare_conical.py`**: Windows file names are case-insensitive, so the DAVID `Unslanted_T` and the BENCH `unslanted_T` write the same `lumerical_conical_unslanted_T.csv`, `rcwa_conical_unslanted_T.fsp` and `figures/conical_unslanted_T.png`. The first DAVID run overwrote the BENCH table; the BENCH case was re-run and reproduced its first numbers exactly (rms 5.4 × 10⁻¹² and 2.1 × 10⁻¹¹). **The DAVID `Unslanted_T` files are now `lumerical_conical_Unslanted_T_DAVID.csv`, `rcwa_conical_Unslanted_T_DAVID.fsp` and `figures/conical_Unslanted_T_DAVID.png`**; the `case` column inside still reads `Unslanted_T`. This is the one deliberate departure from "do not rename files": it is forced by the file system, and the benchmark PC (also Windows) would hit the same collision.
4. Labels in the driver docstring and `AMP_KEYS` comment: conical table mode, the (θ, φ) map, powers and s→p conversion are VERIFIED 2026-10-23; amplitude power normalisation VERIFIED 2026-10-23; the s–p source phase stays UNTESTED.

No fallback was needed for step 2: the RCWA region of the shipped model ("2-D Y-normal") accepted non-zero φ, and `incident angles actual` read back the requested (θ, φ) on every solve (the driver refuses otherwise).

**Side effect to know about.** The Step 0 command `python test_conical.py` itself rewrites `evidence/conical_gates.txt` (line 271 of the script). That file now holds this PC's run of the same 12 gates (0 failures, 362 s, 2026-10-23) instead of the benchmark PC's copy (SHA-256 `7afbee97…ddeb16`). Nothing else under `reference/`, `lumerical/ansys_example/` or `results/reference_conical.csv` was touched (hashes re-checked). Keep the benchmark PC's original when merging.

## Evidence trail

Seven sweeps at 3 harmonics, 160 layers: 189 + 49 + 196 + 98 + 23 + 23 + 23 = 601 tabulated rows (each run also solves five untabulated planar rows for the sign fit). Two driver defects above were found by the first `Slanted_T` and `Unslanted_T` runs, which wrote empty `*.FAILED.csv` files; those were set aside after the fixed re-runs succeeded and are not in the pack. Solver time per sweep: about 11 min (`lumerical_R`), 3 min (`lumerical_T`), 12 min (`zemax_kba_R`), 6 min (`unslanted_T`), 2 min (`Slanted_T`), 2 min (`Slanted_R`), 4 min (`Unslanted_T`).

## Lumerical version and licence

Lumerical 8.35.4609 (Ansys v261). The RCWA solver ran without any licence dialog or error. `ansyslmd.lic` lists the features `lumerical_solve` and `lumerical_gui` (expiry 15 July 2027); the licence edition name is not stated in the licence file or readable through `lumapi`: ASSUMED that it is the standard Lumerical licence.

## Still UNTESTED

- The s–p relative phase of the two sources (the Jones C5 shortfall above) and hence the closed-form corrected model against Lumerical amplitudes.
- `p_source_sign`: +1 chosen by the smaller error on every row, not established independently.
- Layer and harmonics convergence for the six cases other than `lumerical_R` (only the agreement with the referee at 3 harmonics, 160 layers is VERIFIED).
- The 3-D region fallback (never needed).
- `lumerical_T` beyond the Bragg row: grazing signal order; the large closed-form error is expected (outside the domain), the solver-against-referee test passes.
