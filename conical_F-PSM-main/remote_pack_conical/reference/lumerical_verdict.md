# Lumerical RCWA leg of the F-PSM benchmark — verdict

*Written 2026-09-18 on a PC with Ansys Lumerical 2026 R1 (v261, build 4522; FDTD engine reports
`8.35.4522`), driven headless through `lumapi` from
`C:\Program Files\Lumerical\v261\api\python\lumapi.py`, Windows, Python 3.13.14, numpy 2.5.1,
matplotlib 3.11.1. Every number below is VERIFIED 2026-09-18 unless labelled otherwise. British
English. Efficiencies are fractions of 1 W incident, TE ("S"), index-matched, all angles inside
the medium.*

## Verdict

**Lumerical's RCWA solver reproduces the Python RCWA referee to rms 1.0×10⁻¹²–4.3×10⁻⁵ on all four
benchmark gratings and both sweeps of each — gate L1 (≤ 1×10⁻³) passes with two to eight orders of
margin on every one of the eight sweeps, and on the optional fifth case too. It passes against the
*exact* referee, ε = (n₀ + n₁cos)²: against the ε = n₀² + 2n₀·dn·cos form that the OpticStudio
RCWA DLL was found to implement, Lumerical is 5–100× worse (2.1×10⁻⁴–4.9×10⁻³, and 7.6×10⁻⁵
against 1.6×10⁻¹² on the unslanted case). So the two Ansys RCWA solvers do *not* agree with each
other on the permittivity: Lumerical keeps the dn² terms that the DLL drops, and Lumerical is the
one that matches its own documented index equation n = n₀ + Δn·cos(K·r) — which is also the note's
definition and the referee's. On `unslanted_T_rot5`, the one grating both Ansys solvers can
represent, each is an order of magnitude closer to a different referee form on the same points:
Lumerical 3.7×10⁻⁵ exact against 2.1×10⁻⁴ ε-linear, the DLL 1.9×10⁻⁴ exact against 1.0×10⁻⁵
ε-linear. Three further things follow. (1) The referee is independently confirmed: an Ansys
production RCWA lands on it to 4×10⁻⁵ without any tuning, so the DLL–referee gap reported in the
OpticStudio leg was the DLL's, not the referee's. (2) Lumerical solves the exactly unslanted
grating (`unslanted_T`, K_z = 0) that terminates the OpticStudio server, and returns the note's
test vector T-a as 0.586035 against the referee's 0.586035 — agreement to 1.8×10⁻¹³, i.e. the two
codes are solving the identical eigenproblem in double precision. (3) At the `zemax_kba_R`
cut-off where the DLL departed from both referee forms (θ = 17.0°: DLL 0.001978 against 0.044553),
Lumerical returns 0.044550 — it does not misbehave there. F-PSM sits where the referee sits: its
rms distance from Lumerical equals its distance from the referee, two to five orders of magnitude
below Kogelnik's.**

| case | sweep | n | Python RCWA vs Lumerical rms | max | Kogelnik vs Lumerical rms | F-PSM vs Lumerical rms | Python RCWA (ε = n₀²+2n₀dn·cos) vs Lumerical rms | max |
|---|---|---|---|---|---|---|---|---|
| `lumerical_R` | angle −50…30° | 161 | **2.9×10⁻⁵** | 6.0×10⁻⁵ | 1.0×10⁻² | 2.2×10⁻⁴ | 4.2×10⁻³ | 1.0×10⁻² |
| `lumerical_R` | λ 450–650 nm | 101 | **1.8×10⁻⁵** | 6.1×10⁻⁵ | 1.0×10⁻² | 4.1×10⁻⁵ | 2.4×10⁻³ | 1.0×10⁻² |
| `lumerical_T` | angle −20…50° | 141 | **9.9×10⁻⁶** | 2.6×10⁻⁵ | 1.5×10⁻² | 8.3×10⁻³ ¹ | 1.4×10⁻³ | 1.6×10⁻² |
| `lumerical_T` | λ 400–700 nm | 101 | **3.0×10⁻⁵** | 4.4×10⁻⁵ | 8.7×10⁻³ | 6.3×10⁻³ ¹ | 7.3×10⁻⁴ | 1.1×10⁻³ |
| `zemax_kba_R` | angle −70…25° | 191 | **2.3×10⁻⁵** | 7.4×10⁻⁵ | 2.1×10⁻² | 5.1×10⁻⁵ | 4.8×10⁻³ | 1.9×10⁻² |
| `zemax_kba_R` | λ 450–650 nm | 101 | **2.4×10⁻⁵** | 7.6×10⁻⁵ | 1.7×10⁻² | 6.2×10⁻⁵ | 4.9×10⁻³ | 1.9×10⁻² |
| `unslanted_T` | angle 24…44° | 201 | **1.1×10⁻¹²** | 3.8×10⁻¹² | 1.8×10⁻² | 4.6×10⁻⁷ | 7.6×10⁻⁵ | 1.7×10⁻⁴ |
| `unslanted_T` | λ 400–600 nm | 101 | **1.6×10⁻¹²** | 4.0×10⁻¹² | 2.0×10⁻² | 5.3×10⁻⁷ | 8.5×10⁻⁵ | 1.7×10⁻⁴ |
| `unslanted_T_rot5` ² | angle 29…49° | 201 | **3.7×10⁻⁵** | 1.5×10⁻⁴ | 1.8×10⁻² | 3.7×10⁻⁵ | 2.1×10⁻⁴ | 8.4×10⁻⁴ |
| `unslanted_T_rot5` ² | λ 400–600 nm | 101 | **4.3×10⁻⁵** | 1.5×10⁻⁴ | 1.9×10⁻² | 4.3×10⁻⁵ | 2.3×10⁻⁴ | 8.0×10⁻⁴ |

Source: `results/summary.csv` (columns `python_rcwa_vs_lumerical_rcwa_*`,
`kogelnik_vs_lumerical_rcwa_rms`, `fpsm_vs_lumerical_rcwa_rms`,
`python_rcwa_epslin_vs_lumerical_rcwa_*`), figures `figures/compare_<case>.png` (Lumerical is the
pink triangles; the lower panels carry the deviation). The Kogelnik and F-PSM columns are rms over
the points where F-PSM is defined.

¹ `lumerical_T` is deliberately outside F-PSM's stated domain (the m = 2 neighbour order is open,
condition T5 fails), so F-PSM's 6–8×10⁻³ there is the model's own known limitation and not a
Lumerical result. Against the referee F-PSM scores the same 8.3×10⁻³ on this grating.

² `unslanted_T_rot5` is the optional fifth case: T-a with both recording beams turned by 5°, added
to the pack only because the OpticStudio DLL cannot represent K_z = 0. Lumerical does not need it
(it runs `unslanted_T` itself), and it is **not** the note's vector T-a. It was run last, so that
the two Ansys solvers could be compared on one identical grating — see the next section.

## The two Ansys solvers on one grating

`unslanted_T_rot5` is the only grating both Ansys solvers can represent, so it is the clean
head-to-head. Same grating, same sweeps, same referee, scored by the same `compare.py`:

| | vs referee, exact | vs referee, ε = n₀²+2n₀dn·cos | verdict |
|---|---|---|---|
| **Lumerical RCWA**, angle | **3.7×10⁻⁵** | 2.1×10⁻⁴ | 5.5× closer to the **exact** form |
| **Lumerical RCWA**, λ | **4.3×10⁻⁵** | 2.3×10⁻⁴ | 5.4× closer to the **exact** form |
| OpticStudio RCWA DLL, angle | 1.9×10⁻⁴ | **1.0×10⁻⁵** | 19× closer to the **ε-linear** form |
| OpticStudio RCWA DLL, λ | 2.1×10⁻⁴ | **1.2×10⁻⁵** | 18× closer to the **ε-linear** form |

Each solver is an order of magnitude closer to a *different* referee form, on the same grating and
the same points. That is not a numerical-accuracy difference; it is a definitional one, and it
settles the question the OpticStudio leg raised. Note that n₁ = 0.03 here is the mildest grating in
the benchmark — precisely where the OpticStudio leg found the dn² terms small enough for its gate
Z5 to pass — and the separation is still unambiguous in both directions.

`figures/compare_unslanted_T_rot5.png` shows it directly: the two lower panels carry Lumerical's
deviation from the exact referee below the DLL's, while the DLL's deviation from the ε-linear
referee (dashed) drops below both.

## Which ε form does Lumerical follow? — the discriminating points

The off-Bragg rows of EXECUTION_LUMERICAL.md §6 separate the two RCWA forms from each other and
both from Kogelnik. Lumerical lands on the exact column every time; the −14° row of `lumerical_R`
is the sharpest, where the three candidate answers differ by more than a factor of two:

| case | θ | Lumerical | referee, exact | referee, ε-linear | Kogelnik |
|---|---|---|---|---|---|
| `lumerical_R` | 0° (Bragg) | 0.804914 | **0.804947** | 0.805032 | 0.805403 |
| `lumerical_R` | −28° (Bragg) | 0.804914 | **0.804947** | 0.805032 | 0.805403 |
| `lumerical_R` | **−14°** | **0.023933** | **0.023944** | 0.027415 | 0.014014 |
| `lumerical_R` | 7.5° | 0.100138 | **0.100148** | 0.100388 | 0.084730 |
| `lumerical_T` | 0° (Bragg) | 0.978494 | **0.978518** | 0.978292 | 0.987535 |
| `lumerical_T` | 28° (Bragg) | 0.978502 | **0.978518** | 0.978292 | — |
| `lumerical_T` | 5° | 0.013856 | **0.013852** | 0.013689 | 0.025587 |
| `zemax_kba_R` | −45° and 0° (Bragg) | 0.873758 | **0.873784** | 0.873918 | 0.874565 |
| `zemax_kba_R` | −22.5° | 0.012826 | **0.012827** | 0.011815 | 0.042365 |
| `unslanted_T` | 33.8627° (Bragg) | 0.586035 | **0.586035** | 0.585890 | 0.585763 |
| `unslanted_T` | 36° | 0.243709 | **0.243709** | 0.243739 | 0.253146 |

At −14° Lumerical is 1.0×10⁻⁵ from the exact referee, 3.5×10⁻³ from the ε-linear one and
9.9×10⁻³ from Kogelnik. The `unslanted_T` rows are reproduced to 1.8×10⁻¹³, which is
double-precision agreement, not a fit.

### The `zemax_kba_R` cut-off, where the OpticStudio DLL failed

The first order goes evanescent past θ = 17.07°. The OpticStudio leg reported that the DLL
departed from both referee forms at the last two points before that. Lumerical does not:

| θ | Lumerical | referee, exact | referee, ε-linear | OpticStudio DLL |
|---|---|---|---|---|
| 16.0° | 0.018724 | 0.018726 | 0.017974 | 0.017974 |
| 16.5° | 0.014482 | 0.014483 | 0.016592 | 0.016593 |
| **17.0°** | **0.044550** | **0.044553** | 0.049301 | **0.001978** |
| 17.5° | 0.000000 | 0.000000 | 0.000000 | 0.000000 |

Lumerical tracks the exact referee to 3×10⁻⁶ through the cut-off. Note also that 16.0° and 16.5°
are points where the DLL agreed with the ε-linear referee to six decimals — further evidence that
the DLL's expansion, not its cut-off handling alone, is what separates it from Lumerical.

## Converged settings and the acceptance gates

Settings used for every production sweep: **`--harmonics 3 --layers 160`**, i.e. `k vectors
domain` = rectangular with `max number ku` = 3 (orders −3…+3 along x), `max number kv` = 0, and
160 homogeneous layers per grating z-period.

- **L0** (one solve at Bragg) — PASS on all four cases. One order carries the peak
  (`lumerical_R` 0.804914 in order +1, `lumerical_T` 0.978494 in +1, `zemax_kba_R` 0.873758 in +1,
  `unslanted_T` 0.586035 in −1); Σ over all orders of (R + T) = 1.000000 and the solver's own
  `total_energy` R + T = 1.000000 in every case, against a tolerance of 1×10⁻³; the simulated
  thickness reads back as 6.0000 / 6.0000 / 5.0000 / 10.0000 µm as asked.
- **L2** (convergence) — PASS. `results/lumerical_converge_<case>.txt` for all four cases.
  Harmonics are converged at the *first* rung: 3, 5, 9 and 13 agree to all six printed decimals on
  every case at Bragg and off Bragg, so 3 is used. Layers per z-period converge from below:
  for `lumerical_R` at Bragg 10 → 0.796427, 20 → 0.802840, 40 → 0.804422, 80 → 0.804816,
  160 → 0.804914, so the 80 → 160 step is 9.8×10⁻⁵ < 1×10⁻⁴ and 80 is the smallest qualifying
  rung; 160 was used for margin and leaves only 3.3×10⁻⁵ against the referee. `zemax_kba_R` and
  `lumerical_T` behave the same (80 → 160 steps of 7.7×10⁻⁵ and 6.6×10⁻⁵). `unslanted_T` is
  *exactly* converged at every rung (0.586035 at 10 layers as at 160), as it must be: with K_z = 0
  the index does not vary along z at all, so there is no staircase error to converge away.
- **L2, thickness** — PASS. 6 µm is 30.8119 z-periods of 0.194730 µm, and the solver confirms the
  non-integer count: `layer repetition multiple` reads back 30.8119, not 30 or 31. Built the other
  way — rectangle = the whole 6 µm plate, 308 layers, `layer repetition` = "none" — Bragg gives
  0.796426 against the layer-repetition run's 0.796427 at the same 10 layers per period, a
  difference of **1×10⁻⁶**. There is no whole-period truncation.
- **L1** (overlay and score) — PASS against the exact referee on all eight sweeps of the four
  required cases, and on both sweeps of the optional fifth; see the tables above. Expected
  accompanying values are met: `fpsm_vs_lumerical_rcwa_rms` equals F-PSM's distance from the
  referee, and `kogelnik_vs_lumerical_rcwa_rms` is 0.9–2.1×10⁻², in the predicted 1–2×10⁻² band.

The fifth case, `unslanted_T_rot5`, used the same converged settings and passed L0 and L1 the same
way (peak in order −1, energy 1.000000, θ sign +1, rms 3.7×10⁻⁵). No separate convergence ladder
was run for it: it is the same grating as `unslanted_T` at the same n₁ and d, merely rotated, and
its agreement with the referee is already at the level the ladders predict. That ladder is
**UNTESTED** for this case.

Every angle grid the solver returned is exactly `cases.py`'s grid (161 / 141 / 191 / 201 points,
checked to 1×10⁻⁶ degrees), so `compare.py` interpolated nothing. The θ sign was not assumed: the
driver fits both signs and every candidate order against the referee, and chose **+1** on all four
cases, i.e. Lumerical's θ is `cases.py`'s θ directly.

## The NAME MAP: what the probe actually showed

The driver shipped with the pack was written without a Lumerical installation and its guesses were
labelled KNOWN or ASSUMED. `results/lumerical_probe.txt` (regenerated at the end of the run, so it
matches the final driver) lists every property of `::model`, `RCWA` and `rectangle` *with its
current value*, plus every available result and its dataset structure. Against that:

- **`layers` — the property does not exist.** This was one of the two names the brief flagged as
  genuinely unknown. The RCWA object's layer boundaries are its *interfaces*, generated from object
  faces by `interface reference positions`, a table of `[object, face, count]` rows. The shipped
  file holds `[["::model::rectangle","max",10],["::model::rectangle","min",1]]` and reports 11
  evenly spaced `interface actual positions` across the rectangle's z span — 10 layers per
  z-period, exactly the number the article quotes. VERIFIED by experiment: setting the "max" row's
  count to N gives N + 1 positions spanning the same z range (checked at N = 4 → 5 and N = 20 → 21),
  and `layer repetition start/end interface` = 1/2 index the *rows* of that table, so the
  repetition keeps spanning the whole rectangle whatever N is. `set_layers` uses this.
- **The incident-angle sweep** — the other flagged unknown — is `incident angle` = `"range"` with
  `minimum theta` / `maximum theta` / `theta points` / `theta spacing`, or `incident angle` =
  `"table"` with an N × 2 `incident angle table` of (θ, φ). Both are used: range for the uniform
  production sweeps, table for the convergence ladders' three specific angles. `incident angles
  actual` reads back what will really be solved and the driver checks it against what was asked.
  The guesses `theta min` / `theta max` / `number of theta` / `phi` were all wrong.
- **There is no `polarization` property.** The RCWA object solves s and p together and reports
  them as separate result attributes. TE is therefore simply the `Rs_`/`Ts_` attributes; the
  draft's `polarization = "S"` had nothing to set.
- **The result is `grating_power`, not `grating_characterization`.** `grating_characterization`
  carries complex *amplitudes* (`Rss`, `Tss`, …); `grating_orders` carries every retained order
  including evanescent ones; `grating_power` carries the power in each *propagating* order and is
  what the article's figure "Rs power order +1 vs Theta" plots. Layout is
  `[lambda, theta, n, m, forward|backward]` — n the x order, m the y order (a single 0 for this
  2D Y-normal region), and the last axis selected by the dataset's own `forward` = [1, 0] vector.
  The draft's `Rs_grating` / `grating_orders` / 3-D reshape were all wrong.
- Confirmed as shipped: `::model`'s six user properties; the rectangle's `index` and `z span`; the
  RCWA object's `x span` / `z span` / `index`; the whole `layer repetition *` family quoted from
  article 36965778575379; `k vectors domain`, `max number k vectors`, `max number ku`, `max number
  kv`; `wavelength center`; `frequency points`; `simulation z span`; `setup script`.

`NAMES_VERIFIED = True` is set, and no entry of `P` is still ASSUMED.

## Deviations from the procedure, and every change made to the driver

The driver was UNTESTED at the solver, as the brief said. Seven changes were needed. None of them
touches the physics or any reference table; all are recorded in the code with a
`VERIFIED 2026-09-18` comment next to the behaviour that forced them.

1. **`probe()` extended** to print every property *with its value* and every available result with
   its structure, not just the names the draft happened to guess. Without this the two unknown
   names could not have been settled: neither appears in a list of guesses.
2. **The NAME MAP rewritten** as described above, and `NAMES_VERIFIED` set.
3. **`read_thickness` corrected.** `simulation z span` is the span of the *whole* RCWA region after
   repetition, so it includes the background padding the shipped model leaves around the rectangle
   (region z span = 1.5 dz against a rectangle z span of dz). It read 6.097365 µm for a 6 µm
   plate. The grating is `simulation z span − (RCWA z span − rectangle z span)` = 6.000000 µm
   exactly. The padding is background index n₀, which is the benchmark's own no-Fresnel-face
   condition, so it is physically harmless — but taken at face value it would have failed gate L0
   by 97 nm.
4. **`::model`'s setup script re-runs during the solve.** It owns `rectangle z span`, so in
   `--full-structure` mode it silently restored the rectangle to one z-period *after* `configure`
   had resized it to 6 µm — the run returned 0.002236 at Bragg, which is exactly tanh² of one
   0.194730 µm period, a plausible-looking small number and precisely the kind of silent failure
   the brief warns about. Fixed by blanking the script for that mode once it has done its job
   (`set_model` restores it), and — more importantly — `solve()` now re-checks the thickness
   **after** the run rather than only before it, so no future geometry reset can pass unnoticed.
5. **`mesh cells z` is now scaled with the region.** It samples the index over the whole RCWA
   region, so the shipped value (sized for 1.5 z-periods) would alias the sinusoid badly once the
   region becomes the whole plate. It is set from the region's size in z-periods, which reproduces
   the shipped value exactly in the repetition case. (This was investigated as the suspected cause
   of item 4 and turned out not to be it; the scaling is kept because it is right.)
6. **`minimum theta` is rejected if it would exceed the current `maximum theta`.** `unslanted_T`
   asks for 33.86° while the shipped file still holds 30°, so every one of its convergence rungs
   failed with a type error. `set_excitation` now widens the window before narrowing it.
7. **Two ordering bugs in the draft's own logic**, both of which produced wrong numbers before they
   were found:
   - the wavelength sweep reused the *angle* sweep's order **index** against a different order
     axis (`grating_power` only carries orders that propagate somewhere in that solve). For
     `lumerical_R` it reported 7.5×10⁻⁵ at λ = 550 nm instead of 0.804914. The order is now
     matched by its integer value.
   - `choose_sign` compared one fixed integer order across the whole sweep, which is ill-posed
     wherever the K·k order-swap falls *inside* the sweep. On `lumerical_T` (swap at θ = 14°) it
     chose the wrong θ sign and reported rms 2.6×10⁻², failing gate L1 at 3.0×10⁻². It now builds
     the candidate curve with the same swap rule the referee and the output rows use, and selects
     sign +1, order +1, rms 9.9×10⁻⁶.
   In `converge`, `np.allclose(...) is False` was always False (a numpy bool is not the `False`
   singleton), so the interpolation branch was dead; and each ladder inherited the other ladder's
   last setting. Both ladders now set both numbers explicitly.

**Provenance of the tables.** `lumerical_T` and `unslanted_T_rot5` were produced by the final
driver exactly as shipped. `lumerical_R`, `zemax_kba_R` and `unslanted_T` were produced one
revision earlier, before the `choose_sign` refactor of item 7. That refactor is a **provable
no-op** for those three: it only changes the candidate curve where the K·k order swap falls inside
the sweep, and the swap occurs at 0 of 161, 0 of 191 and 0 of 201 points respectively (it needs
cos(θ + ψ) < 0, and their ψ are 14°, 22.5° and −90° against sweeps that never reach it). Where
there is no swap, `first_order_curve` returns `eta_by_order[k]`, which is character-for-character
the expression the earlier revision used, and the row-writing loop was refactored only to share
`conj_index`. `lumerical_T` is the sole case with a swap (72 of 141 points) and is exactly the one
that was re-run. VERIFIED 2026-09-18 by enumerating the swap mask for all five cases.

Deviations from the written procedure itself:

- **The full-structure cross-check was run at 10 layers per z-period (308 layers), not 160.** At
  160 (4 930 layers) a single Bragg point had not returned after 19 minutes and was stopped. Ten
  per period is the figure the brief itself quotes for this check ("308 layers for `lumerical_R`"),
  and the comparison is like-for-like — both sides at 10 — so the gate is met as specified. The
  behaviour of the full-structure build at 160 layers per period is **UNTESTED**.
- **No GUI step was needed anywhere.** `lumapi` accepted blanking `::model`'s setup script, so the
  `unslanted_T` bypass ran entirely through the API and `lumerical/unslanted_T_base.fsp` was never
  created. The grating was *not* rotated: `unslanted_T` is the note's exact T-a vector.
- **No licence obstacle arose.** Sessions opened headless with `hide=True` throughout; no licence
  dialog appeared and the RCWA solver was available for roughly forty solves across the run.
  Licensing is the shared Ansys installation (`ANSYSLIC_DIR` =
  `C:\Program Files\ANSYS Inc\Shared Files\Licensing`). The seat's *edition* string is not exposed
  through `lumapi` and is therefore **not determined** — what is verified is that the RCWA solver
  and the layer-repetition feature are both licensed and functional on this machine.

Nothing in `lumerical/ansys_example/`, `results/zemax_*.csv` or `reference/` was modified.

## What remains UNTESTED

- The full-structure (`--full-structure`) path at more than 10 layers per z-period, as above.
- `--full-structure` on any case other than `lumerical_R`.
- The `"multiple"` mode of `layer repetition` (only `"thickness"` and `"none"` were used).
- P polarisation, oblique φ, and any wavelength-dependent material: all outside this benchmark.
- `use li factorization` was left at the shipped 0. For TE on a 1-D grating Li factorisation is not
  needed (the ladders converge cleanly and the referee is reproduced), but its effect was not
  measured. ASSUMED irrelevant here.
- The ε-form finding is established by ten sweeps on five gratings with n₁ = 0.03–0.05. It has
  not been tested at large n₁, where the dn² terms would be larger still.
- No convergence ladder was run for `unslanted_T_rot5` (see above).

## Step 0, for the record

`python check_env.py` → `READY`. `python test_fpsm_vectors.py` → `0 failure(s)` (81 checks).
`python test_engines.py` → `0 failure(s)` with the expected `SKIP … PSMGRATING_PATH` line.
`python compare.py lumerical_R` reproduced `python_rcwa_epslin_vs_zemax_rcwa_rms` = 2.95×10⁻⁵ on
the angle sweep, so the referee, the scoring and the reference tables are intact on this machine.

## One note for whoever merges this

`RETURN.md`'s sanity list says `summary.csv` should end up with "the eight `..._lumerical_rcwa_...`
columns". This pack's `compare.py` emits **six**, not eight:
`python_rcwa_vs_lumerical_rcwa_rms` / `_max`, `kogelnik_vs_lumerical_rcwa_rms`,
`fpsm_vs_lumerical_rcwa_rms`, `python_rcwa_epslin_vs_lumerical_rcwa_rms` / `_max`. All six are
filled for all ten rows (five cases × two sweeps); nothing is missing or `nan`. Either the main
project's `compare.py` adds two more columns that this pack's copy does not, or "eight" is a slip
in `RETURN.md`. Worth a glance after `merge_back.py` runs the main project's own `compare.py` —
it is a bookkeeping discrepancy, not a gap in the results. VERIFIED 2026-09-18.

The other two sanity checks in `RETURN.md` are met: `figures/compare_lumerical_R.png` shows all
three marker sets (circles = OpticStudio native Kogelnik, squares = OpticStudio RCWA DLL,
triangles = Lumerical RCWA), and the question "exact referee or ε-linear?" is answered above —
**the exact one, on 10 sweeps out of 10.**
