# EXECUTION_LUMERICAL.md — the Lumerical RCWA leg of the F-PSM benchmark

*Written 2026-09-18 for Claude Opus 5 (Claude Code) on a PC with **Ansys Lumerical (Enterprise)**
installed. You have no other context; everything you need is in this folder or in this file. Read
the whole file before running anything. British English throughout. This pack is a slice of a
larger project; you bring results back, you do not rewrite the project.*

## 0. Your job in one paragraph

A closed-form diffraction-efficiency model for volume gratings, **F-PSM**, has been benchmarked
against Kogelnik's theory and against rigorous coupled-wave analysis (RCWA) on four gratings taken
from Ansys' own knowledge-base articles. The "RCWA" in the figures is an independent Python
implementation (the *referee*). Two Ansys solvers were meant to check that referee. The first,
OpticStudio's RCWA hologram DLL, has been run (2026-09-18, see `reference/ansys_rcwa_verdict_opticstudio.md`).
The second, **Lumerical's RCWA solver, has not**: no machine so far had Lumerical. **Your job is to
produce the Lumerical RCWA first-order efficiency curves for the four benchmark gratings, score
them against the referee, F-PSM and Kogelnik, and write a short verdict.** One driver script is
provided, `code/lumerical_rcwa.py`. It is **UNTESTED at the solver** (it was written without
Lumerical); everything else in the pack is proven. Expect to fix property names; §6 tells you how
to know when it is right.

**Rules of engagement.** (1) Every claim you write carries a status and a date: `VERIFIED <date>`,
`ASSUMED`, `UNTESTED`. (2) Never fabricate or tune a number to make curves agree — if Lumerical
disagrees with the referee, that *is* the result; find out why and report it. (3) Fail loudly: no
silent zero, no plausible default; a failed run must not overwrite a good table (the driver writes
`*.FAILED.csv` instead). (4) Never hand-edit `.fsp`, `.zar`, `.zmx` or any binary; change models
only through the application or its API (`lumapi`). (5) Verdict first, evidence after. (6) Do not
modify anything in `lumerical/ansys_example/`, `results/zemax_*.csv`, `results/summary.csv`'s
existing columns, or `reference/`. (7) You do not need OpticStudio for anything here; if it is
installed on this PC, leave it alone.

## 1. Background — the physics (short; the long form is `reference/execution_original.md` §1)

A volume holographic grating (VHG): a slab of thickness *d* with index *n*(**r**) = *n*₀ + *n*₁ cos(**K**·**r**),
recorded by two collimated beams **k₁**, **k₂** inside the medium at wavelength λ_r, so
**K** = (2π*n*₀/λ_r)(**k₁** − **k₂**). Replayed near **k₁** it diffracts into **k₂** (the first order).
Beams from the same side: transmission grating; opposite sides: reflection grating. The quantity
compared is the **first-order diffraction efficiency η** (fraction of incident power, 0–1) versus
replay angle θ (in the medium) at fixed λ, and versus λ at the Bragg angle.

Restrictions of the whole benchmark: TE ("S") polarisation, lossless, the medium outside the slab
has the slab's mean index *n*₀ (no Fresnel faces), all angles inside the medium.

The four gratings, all *n*₀ = 1.5 (`code/cases.py` is the single definition; every engine reads it):

| case | source | *n*₁ | *d* | λ_r = λ | recording beams (in-medium angles) | type | slant Ψ |
|---|---|---|---|---|---|---|---|
| `lumerical_R` | Lumerical article 37118679173267; its figure "Rs power order +1 vs Theta" is the plot to reproduce | 0.04 | 6 µm | 550 nm | 0°, 28° from the opposite side | reflection | 14° |
| `lumerical_T` | same article, "T" variant | 0.04 | 6 µm | 550 nm | 0°, 28° same side | transmission | 76° |
| `zemax_kba_R` | Zemax article 42661816933139, `test_example1` | 0.05 | 5 µm | 550 nm | −45°, 0° opposite side | reflection | 22.5° |
| `unslanted_T` | F-PSM note §6.1, case T-a | 0.03 | 10 µm | 500 nm | +33.8627°, −33.8627° | transmission | −90° (unslanted) |

A fifth entry, `unslanted_T_rot5` (T-a with both beams turned by 5°), exists only because the
OpticStudio DLL cannot represent K_z = 0. Lumerical should not need it (§5 step 4 has the bypass);
run it last and only if time allows — it lets the two Ansys solvers be compared on one grating.

In the article's model the plane of incidence is x–z and the model's inputs are "angle 1",
"angle 2" (in the medium), "holo type" R/T, n0, Dn (the cosine **amplitude**), w_rec. The driver
maps every case to those inputs (`MODEL` in `lumerical_rcwa.py`).

## 2. What is in the pack

```
EXECUTION_LUMERICAL.md          this brief
CLAUDE.md                       points Claude Code at this brief
RETURN.md                       what to bring back and how it is merged into the main project
FPSM_Simple_Implementation.pdf  the F-PSM definition (for reference; nothing to do with it here)
code/cases.py                   THE grating definitions (single source of truth)
code/fpsm.py kogelnik.py rcwa_te.py     the three engines (proven; rcwa_te has an eps_form diagnostic, see §7)
code/test_fpsm_vectors.py test_engines.py   their tests (proven)
code/check_env.py               step 0: Python, packages, lumapi
code/lumerical_rcwa.py          THE DRIVER              ** UNTESTED at the solver **
code/compare.py                 figures + summary; overlays/scores every Ansys table it finds (proven)
lumerical/ansys_example/        the Lumerical article's example as downloaded: vhg_R_6um_model.fsp, licence,
                                the Zemax/Speos companions (not needed)
results/zemax_<case>.csv        OpticStudio native Kogelnik tables (proven; drawn as circles)
results/zemax_rcwa_<case>.csv   OpticStudio RCWA-DLL tables (drawn as squares) -- your curve goes beside them
results/summary.csv             the scores so far (compare.py rewrites it, adding your columns)
figures/compare_<case>.png      the figures so far (compare.py rewrites them)
reference/                      the main project's README, the original brief, the OpticStudio verdict,
                                and the text of the three Ansys articles (the HTML pages return 403 to
                                plain fetches; these were read through the help-centre JSON API)
```

## 3. What is missing

`results/lumerical_rcwa_<case>.csv` for the four cases, `results/lumerical_probe.txt`,
`results/lumerical_converge_<case>.txt`, `lumerical/rcwa_<case>.fsp` (saved by Lumerical), the
regenerated figures and `summary.csv`, and `results/lumerical_verdict.md`.

Table schema (one schema for every Ansys table; `compare.py` depends on the first five columns):
`sweep` (`angle`|`wavelength`), `source_tilt_x_deg` (**= −θ**, degrees), `lambda_um`,
`P_entry_side_W` (reflected first order, per 1 W incident), `P_far_side_W` (transmitted first
order). The driver adds `order`, `harmonics`, `layers_per_period` and the raw values of both
integer orders (`P_*_plus_W`, `P_*_minus_W`); extra columns are ignored by `compare.py`.

## 4. Environment

Windows. Python ≥ 3.10 (3.11–3.13; the main project ran on 3.11.9) with `numpy` and `matplotlib`.
**Use a Windows Python, not WSL**: `lumapi` loads Lumerical's own DLLs. `lumapi.py` ships with
Lumerical under `<install>\api\python\`; the driver finds it by scanning
`C:\Program Files\Lumerical\v*\api\python\lumapi.py` and `C:\Program Files\ANSYS Inc\v*\Lumerical\api\python\lumapi.py`,
else set the environment variable `LUMAPI_PATH` to its full path. Never hard-code an install
path. Run every script from inside `code/`. Lumerical must be ≥ **2025 R1** for the layer-repetition
feature the article's model uses (`reference/ansys_article_36965778575379_layer_repetition.txt`);
on an older version use the driver's `--full-structure` flag (rectangle = the whole plate, 308
layers for `lumerical_R`; slower but the article says the two agree).

The Lumerical session is opened headless (`lumapi.FDTD(filename=..., hide=True)`). If a licence
dialog or an error about the RCWA solver licence appears, stop and report it — do not work around
a licence.

## 5. Procedure

**Step 0 — prove the pack on this machine.** `python check_env.py` → `READY`. `python
test_fpsm_vectors.py` → `0 failure(s)` (81 checks). `python test_engines.py` → `0 failure(s)` (a
`SKIP … PSMGRATING_PATH` line is expected). `python compare.py lumerical_R` must regenerate
`figures/compare_lumerical_R.png` and print a row whose column `python_rcwa_epslin_vs_zemax_rcwa_rms`
is 2.95×10⁻⁵ — that proves the referee, the scoring and the reference tables are intact here.

**Step 1 — probe.** `python lumerical_rcwa.py --probe lumerical_R`. It opens the article's model and
dumps, to `results/lumerical_probe.txt`, every property name of `::model`, `RCWA` and `rectangle`
(three different ways, each tried independently), the model's setup script, the list of available
results and the keys/shapes of the grating-characterisation dataset. **Read that file.** Then edit
the NAME MAP `P` at the top of the driver so that every entry is a real property name (each entry
says whether its current value is KNOWN from the article/.fsp or ASSUMED), and set
`NAMES_VERIFIED = True`. The driver refuses to configure anything until you do. Two names matter
most and are genuinely unknown: the property that sets **how many homogeneous layers one z-period
is sliced into** (the article says 10 in the shipped file), and the **incident-angle sweep**
properties. If the property list cannot be obtained through `lumapi`, open the model in the GUI
once, read the names from the RCWA object's Edit dialog (tabs General / Geometry / Interfaces /
Excitation / Solver / Results) and from the rectangle's Material tab, and note VERIFIED-by-GUI.

**Step 2 — one solve at Bragg.** `python lumerical_rcwa.py --bragg lumerical_R`. Gates L0 (§6). It
prints the efficiency of every order at θ = 0°, λ = 550 nm, the energy sum, the thickness the
solver will really simulate (read back from `simulation z span`) and the referee's value.

**Step 3 — converge.** `python lumerical_rcwa.py --converge lumerical_R` → ladders in harmonics
(3/5/9/13) and layers per z-period (10/20/40/80/160) at Bragg and at the off-Bragg reference
angles (§6, gate L2), written to `results/lumerical_converge_lumerical_R.txt`. Pick the smallest
settings whose next rung moves η by < 1×10⁻⁴ and pass them to the sweeps as `--harmonics N
--layers N`. Also run one Bragg solve with `--full-structure` and confirm it equals the
layer-repetition result to < 1×10⁻⁴ (gate L2, thickness).

**Step 4 — the four sweeps.** `python lumerical_rcwa.py <case> --harmonics N --layers N` for
`lumerical_R`, `zemax_kba_R`, `lumerical_T`, `unslanted_T` → `results/lumerical_rcwa_<case>.csv`,
`lumerical/rcwa_<case>.fsp`. Notes per case:
- `lumerical_R` first: the article reports ≈ 0.80 at θ = 0° and −28°; the referee says 0.804947.
- `lumerical_T`: past θ = 14° the "first order" label swaps to the conjugate integer order
  (OpticStudio's convention, which `compare.py` follows for the Python models). The driver applies
  the same rule and keeps both raw orders in the CSV; check the printout says which integer order
  carries each Bragg peak (θ = 0° and θ = 28°).
- `unslanted_T`: K_z = 0, the model's z-period is infinite and its setup script cannot be used.
  The driver's `configure_unslanted` bypasses the script (see its docstring: stops the setup
  script, sets the rectangle's index equation to n₀ + n₁cos(2πx/dx), rectangle z-span = d,
  repetition off). If `lumapi` refuses to blank the setup script, do it once in the GUI and save as
  `lumerical/unslanted_T_base.fsp`; the driver then uses that file. Do not rotate the grating
  instead — Lumerical, unlike the OpticStudio DLL, has no reason to fail at K_z = 0, and the point
  of this case is the note's exact test vector T-a (referee 0.586035 at Bragg).
- The driver does not assume Lumerical's sign of θ: it fits the solved angle curve to the referee
  under both signs and both candidate orders and reports which it used (rms must be < 0.05 for the
  fit to be accepted; the final score is what `compare.py` says).

**Step 5 — overlay and score.** `python compare.py` (all cases). It picks up the new tables
automatically, adds triangles to `figures/compare_<case>.png` with the deviation from the referee
in the lower panels, and fills the columns `python_rcwa_vs_lumerical_rcwa_rms/_max`,
`kogelnik_vs_lumerical_rcwa_rms`, `fpsm_vs_lumerical_rcwa_rms` in `results/summary.csv`.
**Look at every figure you regenerate before calling it done.**

**Step 6 — write up.** Create `results/lumerical_verdict.md`: verdict first, in bold; a table case ×
sweep with referee-vs-Lumerical rms and max, Kogelnik-vs-Lumerical rms, F-PSM-vs-Lumerical rms,
and — this is the question the OpticStudio leg raised, §7 — whether Lumerical sits on the exact
referee or on the ε-linear one; then every deviation from this procedure, every ASSUMED name you
settled and how (quote the probe), the converged settings, Lumerical's version and licence
edition, and everything that remains UNTESTED. Then follow `RETURN.md`.

## 6. Acceptance gates — how you know it is right

Reference values from the converged Python referee (L = 14), TE, per 1 W incident; "R+1 / T+1" is
the first order on the entry / far face; T0 the undiffracted beam. The last column is the referee
with the OpticStudio DLL's permittivity expansion (§7): a rigorous solver must land on one of the
two RCWA columns, never on Kogelnik's.

| case | θ | λ | referee R+1 / T+1 | T0 | Kogelnik | F-PSM | referee, ε-linear |
|---|---|---|---|---|---|---|---|
| `lumerical_R` | 0° and −28° (Bragg) | 550 nm | R+1 = **0.804947** | 0.1947 | 0.805403 | 0.805021 | 0.805032 |
| `lumerical_R` | −14° | 550 | R+1 = **0.023944** | 0.9760 | 0.014014 | 0.023941 | 0.027415 |
| `lumerical_R` | 7.5° | 550 | R+1 = **0.100148** | 0.8995 | 0.084730 | 0.100133 | 0.100388 |
| `lumerical_T` | 0° (Bragg) | 550 | T+1 = **0.978518** | 0.0139 | 0.987535 | 0.986288 | (see summary.csv) |
| `lumerical_T` | 5° | 550 | T+1 = **0.013852** | 0.9790 | 0.025587 | 0.014074 | |
| `zemax_kba_R` | −45° and 0° (Bragg) | 550 | R+1 = **0.873784** | 0.1256 | 0.874565 | 0.873949 | |
| `zemax_kba_R` | −22.5° | 550 | R+1 = **0.012827** | 0.9868 | 0.042365 | 0.012826 | |
| `unslanted_T` | 33.8627° (Bragg) | 500 | T+1 = **0.586035** | 0.4139 | 0.585763 | 0.586035 | 0.585891 |
| `unslanted_T` | 36° | 500 | T+1 = **0.243709** | 0.7562 | 0.253146 | 0.243707 | |

The off-Bragg rows discriminate: at −14° a genuine RCWA gives 0.0239 (or 0.0274 with the ε-linear
expansion), Kogelnik 0.0140. Empty cells: compute them with
`rcwa_te.signal_efficiency(..., eps_form="linear")` (the `--converge` printout does this for you).

- **L0** (step 2) one order carries R ≈ 0.80 at Bragg for `lumerical_R`; Σ over all orders of
  (R + T) = 1.000 ± 1×10⁻³; the simulated thickness reads back as 6.000 µm. *Else:* geometry or
  name mapping is wrong — check the sign of angle 2, holo type, Dn (amplitude, not peak-to-peak),
  the repetition thickness units (metres).
- **L2** (step 3) harmonics and layers per z-period each converged to < 1×10⁻⁴ at Bragg *and* at the
  off-Bragg angles; layer-repetition thickness equals the full-structure result to < 1×10⁻⁴ (6 µm
  is 30.8 z-periods: a solver that truncates to whole periods would be 0.16 µm short, which moves
  η by ~10⁻³ and shifts the sidelobes).
- **L1** (step 5) `python_rcwa_vs_lumerical_rcwa_rms` ≤ 1×10⁻³ on every case, **or**
  `python_rcwa_epslin_vs_lumerical_rcwa_rms` ≤ 1×10⁻³ if Lumerical turns out to expand ε as the
  OpticStudio DLL does. Report which. Expected then: `fpsm_vs_lumerical_rcwa_rms` ≈ the referee's
  distance, `kogelnik_vs_lumerical_rcwa_rms` ≈ 1–2×10⁻².
- **If a gate cannot be met,** do not loosen it silently: report the measured value, your
  diagnosis, and mark the affected table `UNTESTED` or `FAILED`.

## 7. Pitfalls already paid for (VERIFIED 2026-09-17/18 unless marked)

- **The OpticStudio RCWA DLL implements ε = n₀² + 2n₀·dn·cos(K·r), dropping the dn² terms** of the
  exact ε = (n₀ + dn·cos)², although its article documents n = n₀ + dn·cos. That costs 10⁻³–10⁻² in
  the sidelobes (3.5×10⁻³ at the −14° minimum of `lumerical_R`) and was the whole of the
  DLL–referee gap: with `eps_form="linear"` the referee matched the DLL to rms 2–3×10⁻⁵ on every
  sweep. `compare.py` scores your tables against **both** forms (columns `..._epslin_...`). Which
  one Lumerical follows is a result in itself; do not prejudge it. Lumerical's index equation is
  the exact n = n₀ + Δn·cos, so an exact ε = n² is the natural expectation.
- **Order labelling** (`lumerical_T`, θ > 14°): see step 4. `compare.py` uses the rule
  cos(θ + Ψ) ≥ 0 → the Bragg order, else its conjugate.
- **Cut-off.** Near a first order's cut-off (`zemax_kba_R`, θ = 16.8–17.0°, evanescent past 17.07°)
  the OpticStudio DLL departed from both referee forms (0.002 vs 0.045 at 17.0°). Watch the same
  two points in Lumerical and report them separately if they misbehave; do not drop them silently.
- **Thickness is not in the model**: it is the RCWA object's layer repetition. Read back
  `simulation z span`. `getnamed("RCWA","simulation z span")` is only meaningful while repetition
  is not "none".
- **Dn is the cosine amplitude** (the article: "modulation … ±0.02 … since Δn = 0.04" describes
  the index preview's half-range; the equation is n₀ + Δn·cos, amplitude Δn = 0.04 — the referee
  uses n₁ = 0.04 and reproduces the article's 0.80 peak; keep Dn = 0.04).
- **The angle grid.** RCWA in Lumerical sweeps θ as a range (min/max/number) — the driver sets
  exactly `cases.py`'s grid so `compare.py` needs no interpolation; if the solver returns a
  different grid, `compare.py` interpolates the referee onto it (that is fine, say so).
- **Ansys knowledge-base pages return HTTP 403** to plain fetches; they are readable through
  `https://optics.ansys.com/api/v2/help_center/en-us/articles/<id>.json` with a browser
  user-agent. The three relevant ones are already in `reference/`. Article ids: 37118679173267
  (this model), 36965778575379 (layer repetition), 42661761082899 (dynamic Lumerical–OpticStudio
  link: RCWA object property vocabulary), and the RCWA solver object page linked from them.
- **Long solves**: the article warns computation time "varies greatly" with the angular/spectral
  sampling. The angle sweeps are 141–201 points, the wavelength sweeps 101; if a sweep takes more
  than ~30 min, run it in the background and watch its log rather than waiting.

## 8. Deliverables checklist

- [ ] step 0 outputs: `READY`, `0 failure(s)` twice, `compare.py lumerical_R` reproduces 2.95×10⁻⁵
- [ ] `results/lumerical_probe.txt`; the NAME MAP edited, `NAMES_VERIFIED = True`, every ASSUMED entry either confirmed (say how) or still standing and explained
- [ ] `results/lumerical_converge_<case>.txt` for at least `lumerical_R`; converged settings stated
- [ ] `results/lumerical_rcwa_{lumerical_R,zemax_kba_R,lumerical_T,unslanted_T}.csv` + `lumerical/rcwa_*.fsp`
- [ ] four regenerated `figures/compare_*.png` showing the Lumerical triangles; `results/summary.csv` with the Lumerical columns filled
- [ ] `results/lumerical_verdict.md`
- [ ] the final `code/lumerical_rcwa.py`
- [ ] (optional) `unslanted_T_rot5` as well, so both Ansys solvers meet on one grating

Finish with a short report: verdict; the agreement numbers; which ε form Lumerical follows; what
you changed in the driver and why; what is still untested. If something is blocked by a licence,
an access control or a missing product, say exactly what and stop there — do not work around it.
