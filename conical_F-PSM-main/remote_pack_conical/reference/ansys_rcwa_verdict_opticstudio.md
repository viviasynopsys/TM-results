# Ansys RCWA leg of the F-PSM benchmark — verdict (OpticStudio DLL; Lumerical addendum at the end)

*Written 2026-09-18 on a PC with Ansys Zemax OpticStudio 2024 R1.00 (PremiumEdition seat, standalone
ZOS-API through pythonnet 3.1.0, Python 3.11.9, numpy 2.4.6, matplotlib 3.11.2). Lumerical is not
installed on this PC; the Lumerical leg was run on a second PC on 2026-09-18 and merged here on
2026-09-30 — see the addendum §"Lumerical leg" at the end and `results/lumerical_verdict.md`. Every
number below is VERIFIED 2026-09-18 unless labelled otherwise. British English. Efficiencies are
fractions of 1 W incident, TE, index-matched.*

## Verdict

**OpticStudio's RCWA hologram DLL (`vhg_constant_hologram230313.dll`) runs on this Premium seat and,
once its shipped settings are corrected, reproduces the Python RCWA referee to rms 2–3×10⁻⁵ on all
three benchmark gratings it can represent — but only after one definitional difference is
accounted for: the DLL expands the permittivity as ε = n₀² + 2n₀·dn·cos(K·r), dropping the dn²
terms, although its own article documents n = n₀ + dn·cos(K·r). Against the referee that follows the
documented (and the note's) definition exactly, the DLL deviates by rms 1.4–5.7×10⁻³ and up to
1.6×10⁻² in η on the three n₁ = 0.04–0.05 gratings — so gate Z5 (≤ 1×10⁻³) FAILS as specified
there, and the failure is the DLL's, not the referee's; on the n₁ = 0.03 rotated substitute the
dn² terms are small enough that Z5 is met (rms 1.9–2.1×10⁻⁴). F-PSM sits where the referee sits:
its rms distance from the DLL equals the referee's, 2–3 orders of magnitude below Kogelnik's. The DLL cannot represent the exactly unslanted grating
(`unslanted_T`, K_z = 0): it terminates the OpticStudio server. A 5°-rotated substitute was run and
is reported separately.**

| case | sweep | n | Python RCWA vs DLL rms | max | Kogelnik vs DLL rms | F-PSM vs DLL rms | Python RCWA (ε = n₀²+2n₀dn·cos) vs DLL rms | max |
|---|---|---|---|---|---|---|---|---|
| `lumerical_R` | angle −50…30° | 161 | 4.2×10⁻³ | 1.0×10⁻² | 1.2×10⁻² | 4.2×10⁻³ | **3.0×10⁻⁵** | 6.0×10⁻⁵ |
| `lumerical_R` | λ 450–650 nm | 101 | 2.4×10⁻³ | 1.0×10⁻² | 1.1×10⁻² | 2.4×10⁻³ | **1.8×10⁻⁵** | 6.1×10⁻⁵ |
| `zemax_kba_R` | angle −70…25° | 191 | 5.7×10⁻³ | 4.3×10⁻² | 2.3×10⁻² | 5.9×10⁻³ | 3.4×10⁻³ ¹ | 4.7×10⁻² ¹ |
| `zemax_kba_R` | angle, θ = 16.5°, 17.0° excluded ¹ | 189 | 4.8×10⁻³ | 1.9×10⁻² | — | — | **2.4×10⁻⁵** | 7.3×10⁻⁵ |
| `zemax_kba_R` | λ 450–650 nm | 101 | 4.9×10⁻³ | 1.9×10⁻² | 1.9×10⁻² | 4.9×10⁻³ | **2.4×10⁻⁵** | 8.0×10⁻⁵ |
| `lumerical_T` | angle −20…50° ² | 141 | 1.4×10⁻³ | 1.6×10⁻² | 1.5×10⁻² | 8.1×10⁻³ | **2.5×10⁻⁵** | 5.6×10⁻⁵ |
| `lumerical_T` | λ 400–700 nm | 101 | 8.2×10⁻⁴ | 1.2×10⁻³ | 8.8×10⁻³ | 6.1×10⁻³ | **5.8×10⁻⁵** | 8.4×10⁻⁵ |
| `unslanted_T` | both | — | **FAILED** — the DLL terminates OpticStudio (K_z = 0, infinite z-period; §"Deviations" 4) | | | | | |
| `unslanted_T_rot5` ³ | angle 29…49° | 201 | 1.9×10⁻⁴ | 7.1×10⁻⁴ | 1.8×10⁻² | 1.9×10⁻⁴ | **1.0×10⁻⁵** | 3.7×10⁻⁵ |
| `unslanted_T_rot5` ³ | λ 400–600 nm | 101 | 2.1×10⁻⁴ | 6.8×10⁻⁴ | 1.9×10⁻² | 2.1×10⁻⁴ | **1.2×10⁻⁵** | 3.7×10⁻⁵ |

Source: `results/summary.csv` (columns `python_rcwa_vs_zemax_rcwa_*`, `kogelnik_vs_zemax_rcwa_rms`,
`fpsm_vs_zemax_rcwa_rms`, `python_rcwa_epslin_vs_zemax_rcwa_*`), figures `figures/compare_<case>.png`.
The "Kogelnik vs DLL" and "F-PSM vs DLL" columns are rms over the points where F-PSM is defined
(not past cut-off). Reference: the same models against the exact referee are Kogelnik 1.0–2.1×10⁻²,
F-PSM 2.4×10⁻⁵ (`lumerical_R`), 5.6×10⁻⁵ (`zemax_kba_R`), 8.3×10⁻³ (`lumerical_T`, outside domain),
5.0×10⁻⁷ (`unslanted_T_rot5`).

¹ Two points of `zemax_kba_R` lie within 0.3° of the first order's cut-off (θ = 16.5°, 17.0°; the
order is evanescent past 17.03°). There the DLL returns 0.0166 and 0.0020 where both referee forms
give 0.0166/0.0145 and 0.0493/0.0446 (`results/zemax_rcwa_probe4.txt`); at 15.0–16.0° the DLL still
agrees with the ε-linear referee to 1×10⁻⁶. This grazing-order departure is not the ε-form and is
not explained (UNTESTED further: whether it is the DLL's angular grid meeting the cut-off edge —
a 4001-point grid gave the same 0.0020 — or its z-slicing for a grazing wave). The native Kogelnik
model had a 7×10⁻⁴ discrepancy at the same point on 2026-09-17.

² Past the fringe-parallel angle (θ > 14°) OpticStudio's convention calls the conjugate order "order
1". The DLL's integer orders do not swap; the driver traced both +1 and −1 at every angle and the
table's `P_far_side_W` follows the K·k_playback rule that `compare.py` applies to the Python models.
Both raw orders are kept in the CSV (`P_far_plus_W`, `P_far_minus_W`); e.g. at θ = 20° order +1
carries 0.000234 and order −1 carries 0.034844. Nothing is hidden.

³ `unslanted_T_rot5` is the note's T-a grating with both recording beams turned by 5° (Ψ = −95°,
same |K|, n₁, d; slab = 2.9 z-periods). It is the nearest grating the DLL can represent and is
scored against the referee **on the same rotated geometry**. It is not the note's vector T-a.

## The permittivity finding, and the evidence that it is the DLL's

1. With the shipped settings corrected (below), the DLL converged in z-slices gives, for
   `lumerical_R` at θ = 0°, −14°, 7.5°: **0.805024, 0.027412, 0.100385** (320 layers per period;
   160 layers: 0.804999, 0.027404, 0.100378; 10 layers, the template value: 0.796520, 0.024596,
   0.097872). The exact referee gives 0.804947, 0.023944, 0.100148 (execution.md §6, converged
   L = 14); F-PSM 0.805021, 0.023941, 0.100133; Kogelnik 0.805403, 0.014014, 0.084730.
2. At θ = −14° (the sidelobe minimum) the gap of 3.5×10⁻³ is not thickness (the DLL's thickness
   derivative equals the referee's: 6.0 → 6.036 µm moves η by −2.5×10⁻³ in the DLL and −2.3×10⁻³
   in the referee), not modulation, slant or angle (each would move the Bragg value by ≥ 10⁻³, the
   DLL's Bragg value is within 8×10⁻⁵), not polarisation (the DLL is polarisation-aware: TM gives
   0.727 at Bragg where TE gives 0.797), not harmonics (Max Order 3/5/9/13 identical to 10⁻¹⁰), not
   the angular interpolation grid (0 / 1001 / 4001 / 16001 / −4001 change it by ≤ 3×10⁻⁶).
3. `rcwa_te.solve(..., eps_form="linear")`, which expands ε = n₀² + 2n₀n₁cos(K·r) instead of the
   exact (n₀ + n₁cos)² = n₀² + n₁²/2 + 2n₀n₁cos + (n₁²/2)cos(2K·r), gives **0.805032, 0.027415,
   0.100388** — the DLL's numbers to 1×10⁻⁵ — and over the full sweeps rms 1.8–5.8×10⁻⁵ (table).
   The DC term alone (n₁²/2 in ε, i.e. a mean-index shift of n₁²/(4n₀) = 2.7×10⁻⁴ for n₁ = 0.04)
   accounts for it: keeping the DC term and dropping only the second harmonic changes the exact
   result by ≤ 8×10⁻⁶.
4. The article (42661816933139) documents the DLL's model as n = n₀ + Δn·cos(K·r) (formula image
   `Picture6en-us.jpeg`, read 2026-09-18). The DLL therefore implements a first-order expansion of
   its documented model. It is a small, systematic approximation in the DLL: harmless at the Bragg
   peak (≤ 10⁻⁴), 10⁻³–10⁻² in the sidelobes, largest near sidelobe minima and for thick or
   strongly modulated gratings. Physically it is the difference between "mean index n₀" (the
   note) and "index whose square is the mean permittivity" (the DLL).

Nothing was tuned: the sweeps use dn = n₁ and Holo index = n₀ exactly as the brief specifies; the
ε-linear referee is an added diagnostic column and a thin dashed line, never the reference.

## Deviations from `execution.md`

1. **Lumerical (step 4, gate L1/L2) skipped** — not installed; the owner instructed this. No
   `results/lumerical_rcwa_*.csv`, no `results/lumerical_probe.txt`; `code/lumerical_rcwa.py` is
   untouched and still UNTESTED.
2. **Python version.** The pack asks for ≥ 3.10; the machine's default `python` is 3.8 without the
   packages. Python 3.11.9 was used (`C:\Users\dell\AppData\Local\Programs\Python\Python311`), with
   numpy 2.4.6, matplotlib 3.11.2, pythonnet 3.1.0 installed today. Step 0: 81/81 and all engine
   checks pass (`0 failure(s)` twice, the expected `SKIP … PSMGRATING_PATH` line present); the
   native Kogelnik sweep of `lumerical_R` was re-run on OpticStudio 2024 R1 and agrees with
   `kogelnik.py` to 3.9×10⁻⁸ (angle) and 5.9×10⁻¹³ (wavelength) — gate met.
3. **Sweep settings differ from the template's** (all settled by the probes, `results/zemax_rcwa_probe*.txt`):
   - `Only these orders` 16 → **0**. The shipped value is an order filter: with it, order 0 and every
     transmitted order return zero power (Σ over orders −5…5 = 0.797, gate Z3 fails); with 0 the sum
     is 1.0000 and the order-1 value is unchanged (0.7965196 either way), so no renormalisation.
   - `Layer per period` 10 → **160** (320 for `unslanted_T_rot5`). 10 layers is 8×10⁻³ low at Bragg;
     the error falls as 1/N² (20: 0.80293, 40: 0.80451, 80: 0.80490, 160: 0.80500, 320: 0.80502);
     160 vs 320 ≤ 2.5×10⁻⁵ on `lumerical_R` — gate Z4 met. For the rotated transmission grating
     160 vs 320 is 1.1×10⁻⁴ and 320 vs 640 is 2.8×10⁻⁵ (0.573505 → 0.573477; ε-linear referee 0.573468), so 320 was used there.
   - `Interpolation` 0 kept, **but its meaning is not what the driver assumed.** It is the size of an
     angular lookup grid (article 42661761082899 for the sibling dynamic-link DLL; this DLL's own
     strings say "interp_samp", "Rays' incident angle outside of sampling grids"): 0 or 1 → 201
     points, > 51 → that many, negative → power-conserving interpolation. 16001 changes the result by
     ≤ 3×10⁻⁶ at the probe angles and ≤ 5.3×10⁻⁵ (rms 4×10⁻⁶) over the whole 161-point angle
     sweep (`results/zemax_rcwa_lumerical_R.angle_only_run1.csv`, which is the 16001-grid run), so
     the 201-point grid is converged to the gate. 16001 also takes 36 s per first call and **kills
     the OpticStudio server at every wavelength change** ("Failed to read from an IPC Port: the pipe
     has been ended", `results/zemax_rcwa_probe3.txt`); 401, 2001 and −2001 do not crash.
   - `Max Order` 5 kept (3/5/9/13 identical to 10⁻¹⁰, gate Z4). `Test Mode` 1 kept (0 or 1: no effect
     on any traced value).
4. **`unslanted_T` could not be run on the DLL.** With K_z = 0 the DLL's z-period is infinite and the
   OpticStudio server terminates on the first trace (three separate sessions,
   `results/zemax_rcwa_probe5.txt`); `results/zemax_rcwa_unslanted_T.FAILED.csv` is the empty
   table. Rotating the grating by 0.01°, 0.1°, 1° (z-period 1714, 171, 17 µm — the 10 µm slab is
   thinner than one period) runs but returns 0.884, 0.485, 0.049 at Bragg against 0.586: the DLL
   does not handle a slab thinner than one z-period either. From 2° (1.2 periods in the slab) the
   DLL agrees with the referee on the rotated geometry (2°: 0.584034 vs 0.584035; 5°: 0.573614 vs
   0.573615 exact / 0.573468 ε-linear at 160 layers; `results/zemax_rcwa_probe6.txt`). The 5°
   rotation was added to `cases.py` as `unslanted_T_rot5` and swept; see the table.
5. **Driver changes** (`code/zemax_rcwa_nsc.py`, all labelled in the file): CLI flags `--interp`,
   `--layers`, `--max-order`; the `DLL_OVERRIDES` dictionary (only `Only these orders` = 0); the
   K·k_playback order rule with both integer orders recorded; settings recorded in extra CSV columns
   (`order`, `interpolation`, `layer_per_period`, `max_order`, `P_*_plus_W`, `P_*_minus_W` — the
   schema `compare.py` reads is unchanged); the probe extended with Layer-per-period and
   Interpolation convergence; the DLL log tail printed. Five probe scripts were added
   (`zemax_rcwa_probe2.py` … `probe6.py`) and the optional `zemax_epe_rcwa.py`.
   `code/rcwa_te.py` gained the optional `eps_form` argument (default unchanged; `test_engines.py`
   still 0 failures). `code/compare.py` scores and draws the ε-linear referee as a diagnostic.
6. **Bench files.** `zemax/rcwa_dll_<case>.zmx` (+ `.ZDA`) were written by OpticStudio 2024 R1 for
   the three completed cases and the rotated substitute; `zemax/kogelnik_native_lumerical_R.zmx`
   was rewritten by the step-0 re-run. Nothing in `zemax/ansys_kba/` was modified; the EPE archive
   was restored by `Tools.OpenRestoreZAR` into `zemax/epe_restored/` (outside `ansys_kba/`).
7. **The MCP servers** in the two ZIPs were not used: the pack's own ZOS-API drivers (pythonnet,
   standalone connection) did everything; nothing needed a GUI.

## ASSUMED items in the brief, and how each was settled

| item | settled |
|---|---|
| name ↔ value pairing of the 24 DLL parameters (execution.md §5 step 2) | VERIFIED: the live printout (`results/zemax_rcwa_probe.txt`) lists exactly the 24 names in the stated order with the stated template values; indices 0-based |
| `Holo Thic (um)` in micrometres | VERIFIED at gate Z2: 6 gives the 6 µm response (0.797 at 10 layers → 0.805 converged); 5.84 / 6.04 / 6.2 µm move η like the referee's thickness derivative |
| "+1" is the first order | VERIFIED by measurement at Bragg for every case (the −1 trace carries ≤ 10⁻⁴); for `lumerical_T` the conjugate labelling past 14° is handled as in note ² |
| Interpolation 0 = exact per-ray solve | **WRONG**: it is a 201-point angular grid (see Deviations 3); adequate to 5×10⁻⁵ |
| the template's other settings can be kept as shipped | **WRONG for `Only these orders`** (an order filter, see Deviations 3); right for Test Mode, Stochastic mode, Rotate Grating, Use Coat File |
| the article's index model | VERIFIED from the formula image: n = n₀ + Δn·cos(K·r); the DLL implements ε = n₀² + 2n₀Δn·cos |

## What remains UNTESTED

- ~~The Lumerical leg entirely~~ — run 2026-09-18 on a second PC, merged 2026-09-30; see the addendum.
- Why the DLL departs from both referees within 0.3° of first-order cut-off (note ¹).
- Why `Interpolation` ≥ 16001 kills the server (a memory or grid-allocation failure is a guess).
- Whether the DLL's ε-expansion is deliberate; Ansys has not been asked. The evidence is numerical
  (a 1×10⁻⁵ match over 800+ points with one change to the referee), not a reading of the DLL's code.
- The DLL's behaviour for TM and for n₁/n₂ ≠ Holo index (Fresnel faces) — outside this benchmark.
- The EPE picture is a single trace with the archive's own ray count; no convergence in rays.

## Step 7 (optional) — the article's EPE picture

`zemax/ansys_kba/EPE with RCWA hologram.zar` was restored and traced as shipped (16 objects; object
2 is a Diffraction Grating split by `vhg_constant_hologram230313.dll`, orders −1…1). Detector 14
(300 × 300 px, 6.45×10⁻⁶ W) shows the QR-code slide through the eye lens — the article's headline
picture — saved as `figures/epe_rcwa.png`; detector 16 (300 × 150 px, 1.05×10⁻⁴ W, the eyebox) as
`figures/epe_rcwa_det16_eyebox.png`. **F-PSM was not added and cannot be**: the waveguide index
differs from the hologram index and F-PSM is index-matched only (note §1; README).

## Addendum, 2026-09-30 — the Lumerical leg, and what it settles

The pack built from this project (`../17-09 lumerical pack`, brief `EXECUTION_LUMERICAL.md`) was
run by Claude Opus 5 on a PC with Ansys Lumerical 2026 R1 (v261, FDTD engine 8.35.4522) on
2026-09-18; its tables, probe and convergence printouts and verdict (`results/lumerical_verdict.md`)
were merged into this project on 2026-09-30 and rescored here by the main referee (identical to the
second PC's scores to the third digit on every column).

**Lumerical's RCWA reproduces the exact referee, ε = (n₀ + n₁cos)², on all four benchmark
gratings and both sweeps (rms 1×10⁻¹²–4×10⁻⁵; gate L1 met with two to eight orders of margin),
and the exactly unslanted vector T-a to 1.8×10⁻¹³. Against the ε = n₀² + 2n₀dn·cos form that the
OpticStudio DLL implements, Lumerical is 5–100× worse. So the two Ansys RCWA solvers disagree
with each other on the permittivity; Lumerical is the one that follows the documented index
equation n = n₀ + Δn·cos(K·r), which is also the note's and the referee's; the DLL–referee gap
reported above was the DLL's.** VERIFIED 2026-09-18 (second PC) and 2026-09-30 (rescored here).

| case | sweep | n | Python RCWA vs Lumerical rms | max | Kogelnik vs Lumerical rms | F-PSM vs Lumerical rms | Python RCWA (ε = n₀²+2n₀dn·cos) vs Lumerical rms |
|---|---|---|---|---|---|---|---|
| `lumerical_R` | angle | 161 | **2.9×10⁻⁵** | 6.0×10⁻⁵ | 1.0×10⁻² | 2.2×10⁻⁴ | 4.2×10⁻³ |
| `lumerical_R` | wavelength | 101 | **1.8×10⁻⁵** | 6.1×10⁻⁵ | 1.0×10⁻² | 4.1×10⁻⁵ | 2.4×10⁻³ |
| `zemax_kba_R` | angle | 191 | **2.3×10⁻⁵** | 7.4×10⁻⁵ | 2.1×10⁻² | 5.1×10⁻⁵ | 4.8×10⁻³ |
| `zemax_kba_R` | wavelength | 101 | **2.4×10⁻⁵** | 7.6×10⁻⁵ | 1.7×10⁻² | 6.2×10⁻⁵ | 4.9×10⁻³ |
| `lumerical_T` | angle | 141 | **9.9×10⁻⁶** | 2.6×10⁻⁵ | 1.5×10⁻² | 8.3×10⁻³ | 1.4×10⁻³ |
| `lumerical_T` | wavelength | 101 | **3.0×10⁻⁵** | 4.4×10⁻⁵ | 8.7×10⁻³ | 6.3×10⁻³ | 7.3×10⁻⁴ |
| `unslanted_T` | angle | 201 | **1.1×10⁻¹²** | 3.9×10⁻¹² | 1.8×10⁻² | 4.6×10⁻⁷ | 7.6×10⁻⁵ |
| `unslanted_T` | wavelength | 101 | **1.6×10⁻¹²** | 4.1×10⁻¹² | 2.0×10⁻² | 5.3×10⁻⁷ | 8.5×10⁻⁵ |
| `unslanted_T_rot5` | angle | 201 | 3.7×10⁻⁵ | 1.5×10⁻⁴ | 1.8×10⁻² | 3.7×10⁻⁵ | 2.1×10⁻⁴ |
| `unslanted_T_rot5` | wavelength | 101 | 4.3×10⁻⁵ | 1.5×10⁻⁴ | 1.9×10⁻² | 4.3×10⁻⁵ | 2.3×10⁻⁴ |

At the discriminating point θ = −14° of `lumerical_R`, Lumerical gives 0.023933 (referee
0.023944, ε-linear 0.027415, DLL converged 0.027412, Kogelnik 0.014014). At the `zemax_kba_R`
cut-off point θ = 17.0° where the DLL returned 0.001978, Lumerical returns 0.044550 against the
referee's 0.044553 — the cut-off departure is the DLL's alone. On `unslanted_T_rot5`, the one
grating both solvers can represent, Lumerical is 5.5× closer to the exact form and the DLL 19×
closer to the ε-linear form, on identical points.

Lumerical's converged settings: 3 harmonics (3 = 5 = 9 = 13 to six decimals), 160 layers per
z-period (80 → 160 ≤ 9.8×10⁻⁵), layer-repetition thickness read back as 30.81 periods = 6.000 µm
and confirmed against a full 308-layer structure to 1×10⁻⁶. The unslanted grating is a single
z-homogeneous layer and is converged at any layer count. No convergence ladder was run for
`unslanted_T_rot5` (its Bragg value is 1.4×10⁻⁴ high; UNTESTED further). Six of the pack driver's
guessed property names were wrong and were settled from the solver's own property list
(`results/lumerical_probe.txt`); the final driver is `code/lumerical_rcwa.py`. Nothing was tuned.

What remains UNTESTED after both legs: why the DLL departs near cut-off; why it crashes at
K_z = 0 and at Interpolation ≥ 16001; whether Ansys intends the DLL's expansion (not asked);
large-n₁ and TM behaviour of both solvers.
