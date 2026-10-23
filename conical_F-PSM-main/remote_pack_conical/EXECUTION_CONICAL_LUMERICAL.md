# EXECUTION_CONICAL_LUMERICAL.md — Lumerical RCWA at conical incidence, seven gratings

*For the person at the keyboard: extract the pack into a folder of its own, open Claude Code in the
extracted `remote_pack_conical` folder and type: **"Follow EXECUTION_CONICAL_LUMERICAL.md exactly, step
by step, and give me the report it asks for at the end."** About three hours of solver time.
Written 2026-10-23. British English.*

---

## What you are doing

Seven volume holographic gratings are to be solved by Lumerical's RCWA at **conical incidence**: the
incident wave has a polar angle θ from the plate normal and an azimuth φ about it, so it lies out of
the plane that contains the grating vector. In that geometry s and p polarisations couple: an
s-polarised input diffracts partly into p. You are measuring what Lumerical returns and scoring it
against an independent conical RCWA referee whose values are already in
`results/reference_conical.csv`. You are measuring, not modelling: never tune a number, never invent
one, never continue past a failed gate as if it had passed. Every sentence you write carries
`VERIFIED <date>`, `ASSUMED` or `UNTESTED`.

The cases are defined once, in `code/cases_conical.py`:

| set | case | what | azimuths φ (deg) | rows |
|---|---|---|---|---|
| BENCH | `lumerical_R` | the Lumerical article's reflection grating (n₀ 1.5, Δn 0.04, 6 µm, 550 nm) | 100, 120, 150, 180 | 189 |
| BENCH | `lumerical_T` | its transmission variant | 100 | 49 |
| BENCH | `zemax_kba_R` | the OpticStudio article's reflection grating (Δn 0.05, 5 µm) | 100, 120, 150, 180 | 196 |
| BENCH | `unslanted_T` | the unslanted transmission grating (Δn 0.03, 10 µm, 500 nm) | 30, 45 | 98 |
| DAVID | `Slanted_T` | the author's Table 1 slanted transmission family (Λ 0.342 µm, Ψ 80°, n₁ 0.024, 6.667 µm, 532 nm) | 150 | 23 |
| DAVID | `Slanted_R` | his slanted reflection family (Λ 0.184 µm, Ψ 30°) | 150 | 23 |
| DAVID | `Unslanted_T` | his unslanted transmission family (Λ 0.485 µm, Ψ 90°) | 150 | 23 |

Each (case, φ) is an angle sweep about the conical Bragg angle at the design wavelength. The
BENCH gratings are the ones the September pack already builds in the Lumerical model; the DAVID
gratings are built from two recording beams the driver computes (printed at the start of the run).
Efficiencies are fractions of the incident power; angles are inside the medium (n₀ = 1.5 on both
sides, no Fresnel faces). "s" and "p" are the incident polarisations; each diffracted order's
power sums both output polarisations.

## Before you start

1. Windows (not WSL), Python 3.10 or later with `numpy`, `scipy`, `matplotlib`; Ansys Lumerical
   2025 R1 or later with the RCWA solver licensed; `lumapi.py` under `<Lumerical install>\api\python\`.
   If the driver cannot find it, set `LUMAPI_PATH` to its full path. Never hard-code a path.
2. Work from inside `code\`. Every command below is run there.
3. Do not modify `lumerical\ansys_example\`, `reference\`, `evidence\` or `results\reference_conical.csv`.
4. Do not hand-edit `.fsp` files; the driver changes the model through `lumapi` only.
5. A sweep of ~200 rows at 160 layers per period takes tens of minutes; run it in the background and watch the log.
6. If a licence dialog or an RCWA-licence error appears: stop, report it word for word, do not work around it.

## Step 0 — prove the pack on this machine (10 minutes)

| command | must print | if not |
|---|---|---|
| `python check_env.py lumerical` | `READY` | fix what it names (usually `LUMAPI_PATH`) |
| `python test_conical.py` | `0 failure(s) of 12 gates` (about 90 s) | stop; the pack is damaged |
| `python test_fpsm_vectors.py` | `0 failure(s)` | stop; the pack is damaged |
| `python compare_conical.py lumerical_R` | `no Lumerical table yet` | stop and report |

## Step 1 — probe (5 minutes)

`python lumerical_rcwa_conical.py --probe`

Writes `results\lumerical_probe.txt`. Confirm in it: the RCWA object lists `incident angle`,
`incident angle table`, `incident angles actual`, `minimum phi`, `maximum phi`, `phi points`
(VERIFIED 2026-09-18 on Lumerical v261, FDTD 8.35.4522); result `grating_power` has `Rs_power`,
`Ts_power`, `Rp_power`, `Tp_power`; result `grating_characterization` has `Rss`, `Rsp`, `Rps`, `Rpp`,
`Tss`, `Tsp`, `Tps`, `Tpp`. If a name differs on this build, correct it in the NAME MAP `P` near the
top of `lumerical_rcwa.py` (or `AMP_KEYS` in `lumerical_rcwa_conical.py`) and write down what you
changed and why.

## Step 2 — one conical Bragg point per azimuth (10 minutes)

`python lumerical_rcwa_conical.py --bragg lumerical_R --harmonics 3 --layers 160`

The driver first solves five planar rows (φ = 0 in Lumerical's axes, θ from −3° to +3° about the
planar Bragg angle 0°) and fits the sign of Lumerical's θ and the signal order against the referee
(gate C1 below). Then one solve per azimuth at the conical Bragg angle. Gates:

- **C1** planar fit: `theta sign +1, signal order +1` with `rms vs referee` ≤ 10⁻³ (expected ≈ 3 × 10⁻⁵,
  the September value). A sign of −1 is acceptable if the rms is as small; report it.
- **C0** power: every printed power sum = 1.000 ± 0.001 for s and for p.
- **C3, the conical signature**: at φ = 150°, θ = 24.369°, the s-incident efficiency of the signal
  order is **0.7872** and almost all of it is p-polarised on exit (s→p **0.7872**, p→s **0.7516**;
  the table below). A solver that ignored φ would return the planar value 0.805 with no conversion.

If the solver rejects a non-zero φ (an error on `incident angle table`, or `incident angles actual`
reading back zeros in the second column): the RCWA region of the shipped model is "2D Y-normal".
UNTESTED fallback: set the RCWA object's dimension to 3D with `y span` equal to its `x span`, keep
`max number kv` = 0, and rerun step 2. Record exactly what you changed. If that also fails, stop and report.

## Step 3 — convergence (15 minutes)

Run step 2's command three more times with `--harmonics 5 --layers 160`, `--harmonics 3 --layers 80`
and `--harmonics 3 --layers 320`, and tabulate the four Bragg-point values (s and p) per azimuth.
Gate **C4**: harmonics 3 → 5 changes nothing above 10⁻⁵; layers 80 → 160 → 320 converges from below
as 1/N² (differences falling by about 4). Keep 3 harmonics and 160 layers unless the ladder says
otherwise, and say so if it does.

## Step 4 — the seven sweeps (about three hours)

```
python lumerical_rcwa_conical.py lumerical_R  --harmonics 3 --layers 160
python lumerical_rcwa_conical.py lumerical_T  --harmonics 3 --layers 160
python lumerical_rcwa_conical.py zemax_kba_R  --harmonics 3 --layers 160
python lumerical_rcwa_conical.py unslanted_T  --harmonics 3 --layers 160
python lumerical_rcwa_conical.py Slanted_T    --harmonics 3 --layers 160
python lumerical_rcwa_conical.py Slanted_R    --harmonics 3 --layers 160
python lumerical_rcwa_conical.py Unslanted_T  --harmonics 3 --layers 160
```

Each writes `results\lumerical_conical_<case>.csv` and `lumerical\rcwa_conical_<case>.fsp`, and
prints per azimuth `s rms vs referee`, `p rms`, and the power-sum range. A failed run writes
`*.FAILED.csv` and never overwrites a good table; report it, do not retry blindly.

- `unslanted_T` uses the K_z = 0 bypass of the September driver (index equation written directly,
  repetition off); do not rotate the grating instead.
- `lumerical_T` has only one conical Bragg cone within reach, at θ_B = 69.7° (φ = 100°), where the
  signal order is close to grazing (γ₁/β = 0.22 at Bragg, 0.06 at 75°). The report certifies
  nothing for grazing orders and the closed forms are outside their domain there: corrected-model
  error 1.9 × 10⁻³ at Bragg, up to 0.57 at the sweep's far edge (`reference_conical.csv`). That row
  is a solver-against-referee test only; a large closed-form error on it is expected and is not a
  solver fault. Both neighbouring orders are closed on the whole sweep.
- The DAVID gratings are built from recording beams printed at the start of the run; the driver
  checks they reproduce the grating vector to 10⁻⁹ and refuses otherwise.

## Step 5 — score (5 minutes)

`python compare_conical.py`

Prints one line per (case, φ) and writes `results\summary_conical.csv` and `figures\conical_<case>.png`.
Look at every figure before calling it done. Gates on the summary:

- **C2** `rms_s` and `rms_p` ≤ 10⁻³ on every (case, φ); expected 10⁻⁵ to 4 × 10⁻⁵ on the slanted
  gratings, ≈ 10⁻¹² on `unslanted_T` (the planar values of September).
- **C0** `power_*_min/max` within 1 ± 10⁻³.
- **C5** `amplitude_convention` resolved (not "unresolved") on every row; then
  `jones_err_lumerical_vs_referee` ≤ 10⁻³ and `jones_err_corrected_vs_lumerical` within a factor two
  of `corr_max_err_vs_referee` on the same row. If it stays unresolved, the Jones columns are empty:
  report that and the first ten rows of `grating_characterization` for one case so the convention can
  be established on the benchmark PC. Do not guess a normalisation.
- `p_source_sign` is reported as ASSUMED: the sign of Lumerical's p source relative to the referee's
  basis cannot be fixed from powers; the comparison tries both and keeps the smaller error.

## The numbers that decide whether it worked

Conical referee (21 orders), per 1 W incident, at the conical Bragg angle of each (case, φ). The
closed-form columns are the polarisation-uniform error of the report's Eq. (65) against the referee
and are what the Ansys data will be compared with. Full table: `results\reference_conical.csv`.

| case | φ (°) | θ_B (°) | η(s in) | η(p in) | s→p | p→s | corrected model err | leading err |
|---|---|---|---|---|---|---|---|---|
| lumerical_R | 100 | 4.958 | 0.738795 | 0.802714 | 0.738789 | 0.802453 | 7.5e-5 | 5.3e-4 |
| lumerical_R | 120 | 14.212 | 0.753102 | 0.787188 | 0.753064 | 0.785556 | 1.1e-5 | 4.7e-4 |
| lumerical_R | 150 | 24.369 | 0.787246 | 0.753167 | 0.787211 | 0.751631 | 5.7e-6 | 4.5e-4 |
| lumerical_R | 180 | 28.000 | 0.804947 | 0.736870 | 0 | 0 | 7.3e-5 | 4.6e-4 |
| lumerical_T | 100 | 69.712 | 0.831961 | 0.957615 | 0.015403 | 0.023845 | 1.9e-3 | 2.6e-2 |
| zemax_kba_R | 100 | 8.228 | 0.698913 | 0.867735 | 0.698863 | 0.865741 | 2.5e-5 | 9.3e-4 |
| zemax_kba_R | 120 | 23.402 | 0.732685 | 0.823637 | 0.732420 | 0.811703 | 2.7e-5 | 1.4e-3 |
| zemax_kba_R | 150 | 39.468 | 0.824396 | 0.733723 | 0.824165 | 0.723542 | 1.3e-5 | 1.4e-3 |
| zemax_kba_R | 180 | 45.000 | 0.873784 | 0.694665 | 0 | 0 | 1.7e-4 | 7.8e-4 |
| unslanted_T | 30 | 40.046 | 0.486201 | 0.555040 | 0.474075 | 0.474075 | 4.7e-6 | 7.1e-3 |
| unslanted_T | 45 | 52.000 | 0.611891 | 0.235696 | 0.198089 | 0.198089 | 4.6e-6 | 1.0e-2 |
| Slanted_T | 150 | 25.000 | 0.651751 | 0.456301 | 0.294495 | 0.413785 | 1.4e-6 | 9.9e-4 |
| Slanted_R | 150 | 20.000 | 0.579056 | 0.648305 | 0.211835 | 0.276597 | 8.4e-7 | 2.0e-4 |
| Unslanted_T | 150 | 25.000 | 0.668372 | 0.555199 | 0.496373 | 0.496373 | 8.9e-7 | 4.3e-4 |

Off-Bragg rows discriminate as well: at θ_B + 2° on `lumerical_R`, φ = 150°, the referee gives
η(s) = 0.758304, η(p) = 0.700280; on `unslanted_T`, φ = 30°, 0.141704 and 0.084094.

The φ = 180° rows of the reflection gratings are planar replays at the second Bragg angle (the
September tables' θ = −28° and −45°) and must reproduce the September TE/TM values above to the
September accuracy — that is the continuity check between the two packs.

## Step 6 — write the verdict and return

Write `results\lumerical_conical_verdict.md`: the verdict first, in bold; the C0–C5 numbers per
(case, φ) from `summary_conical.csv`; the convergence ladder; the amplitude convention found and
the evidence; every name or assumption you changed; Lumerical version and licence edition; what
remains UNTESTED. Then follow `RETURN_CONICAL.md`.

Finish with a short report in the chat: verdict; the C2 and C5 numbers; what you changed in the
driver and why; what is still untested.

## Pitfalls already paid for

- There is no polarisation property; s and p are always both solved. s = `Rs_power`/`Ts_power`,
  p = `Rp_power`/`Tp_power`.
- Thickness is the RCWA object's layer repetition, not a model property; the setup script re-runs at
  every solve and restores geometry it owns; the driver reads the real thickness back and refuses on
  a mismatch.
- Layers per period = the `interface reference positions` count (shipped 10; use 160).
- `Dn` is the cosine amplitude, not peak-to-peak.
- Lumerical's θ sign is fitted by the driver, never assumed (+1 on all cases in September); the same
  map is applied to φ (sign −1 means φ → 180° − φ). Efficiencies do not depend on the sign of φ.
- `incident angle table` is an N × 2 array of (θ, φ) in degrees; the driver reads
  `incident angles actual` back and refuses on any mismatch.
- The probe dump is printed in full (numpy threshold raised).
- Ansys knowledge-base pages return HTTP 403 to plain fetches; the relevant articles are in `reference\`.
