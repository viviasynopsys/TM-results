# execution.md — close the Ansys RCWA leg of the F-PSM benchmark

*Written 2026-09-17 for an agent with no prior context, working on a PC that has **Ansys Zemax
OpticStudio (Enterprise)** and **Ansys Lumerical** installed. Everything you need is in this folder
or in this file. Read the whole file before running anything. British English throughout.*

## 0. Your job in one paragraph

A closed-form diffraction-efficiency model for volume gratings, **F-PSM**, has been implemented and
benchmarked against Kogelnik's theory and against rigorous coupled-wave analysis (RCWA). The Kogelnik
leg is already tied to Ansys (OpticStudio's native model, agreement 4×10⁻⁸). The **RCWA leg is not**:
the machine that produced this pack could not run either Ansys RCWA solver, so the "RCWA" in every
figure is an independent Python implementation. **You must produce the missing Ansys RCWA curves —
OpticStudio's RCWA hologram DLL and Lumerical's RCWA solver — for the four benchmark gratings,
overlay them on the existing figures, score all models against them, and write the verdict.**
Two driver scripts are provided for this. Both are **UNTESTED** where they touch the Ansys RCWA
solvers (they could not be run), and proven everywhere else. Expect to debug them; §6 tells you how
to know when they are right.

**Rules of engagement.** (1) Every claim you write carries a status and a date: `VERIFIED <date>`,
`ASSUMED`, `UNTESTED`. (2) Never fabricate or tune a number to make curves agree — if an Ansys curve
disagrees with the Python referee, that *is* the result; find out why and report it. (3) Fail
loudly: no silent zero, no plausible default; a failed run must not overwrite a good table (the
drivers write `*.FAILED.csv` instead). (4) Never hand-edit `.zmx`, `.zar`, `.fsp` or any lens/model
binary; change them only through the application or its API. (5) Verdict first, evidence after.
(6) Do not delete or modify anything in `zemax/ansys_kba/` or `lumerical/ansys_example/`.

## 1. Background — the physics

**The object.** A volume holographic grating (VHG): a flat slab of thickness *d* whose refractive
index varies sinusoidally, *n*(**r**) = *n*₀ + *n*₁ cos(**K**·**r** + φ). It is recorded by the
interference of two collimated beams of wavelength λ_r inside the medium with unit wave-vectors
**k₁**, **k₂**; the grating vector is **K** = (2π*n*₀/λ_r)(**k₁** − **k₂**). Replayed near **k₁** it
diffracts into **k₂** (the "first order"). Beams from the same side give a *transmission* grating
(first order leaves the far face); from opposite sides, a *reflection* grating (first order returns
through the entry face). The quantity of interest is the **first-order diffraction efficiency η** —
the fraction of incident power in that order, 0–1 — as a function of replay angle θ and wavelength λ.

**The three models compared.**

| Model | Nature | Known weakness |
|---|---|---|
| **Kogelnik (1969)** two-wave coupled-wave theory | closed form; what OpticStudio's native volume hologram uses | drops second derivatives and all orders but 0 and 1; error ~10⁻² in η even for ordinary gratings |
| **RCWA** | rigorous: Floquet expansion in 2L+1 orders, eigen-solve, boundary matching | slow; must be converged in L (and, in Ansys tools, in layers per period) |
| **F-PSM** | closed form, a two-wave theory corrected perturbatively for the two nearest neighbour orders (m = −1, 2) and for surface terms; defined in `FPSM_Simple_Implementation.pdf` | lossless, TE only, **index-matched exterior only**; has a stated working domain (conditions T1–T6 / R1–R6 in the note) |

**Common restrictions of this whole benchmark:** TE (σ, "S") polarisation — electric field
perpendicular to the plane of incidence; lossless (real *n*₁); the medium outside the slab has the
slab's mean index *n*₀ (no Fresnel reflections at the faces); all angles are measured **inside the
medium**.

**Conventions (get these wrong and nothing agrees).** Axes in the plane of incidence: *n* along the
surface normal pointing into the plate, *t* tangential. Replay wave **k** = β(sin θ, cos θ),
β = 2π*n*₀/λ. F-PSM writes **K** = *K*(−sin Ψ, cos Ψ) with slant Ψ and takes θ_c = θ + Ψ (angle
between replay beam and **K**). `code/cases.py` is the single definition of every grating and
derives all of this from (**k₁**, **k₂**); every engine reads it. In OpticStudio the plane of
incidence is y–z, a Source Ray tilted τ about X points along (0, −sin τ, cos τ), hence **τ = −θ**.

**The four benchmark gratings** (all *n*₀ = 1.5):

| case | source | *n*₁ | *d* | λ_r = λ | **k₁**, **k₂** (in-medium angles) | type | Ψ |
|---|---|---|---|---|---|---|---|
| `lumerical_R` | Lumerical article 37118679173267; its figure "Rs power order +1 vs Theta" is the plot to reproduce | 0.04 | 6 µm | 550 nm | 0°, 28° from the opposite side | reflection | 14° |
| `lumerical_T` | same article, "T" variant | 0.04 | 6 µm | 550 nm | 0°, 28° same side | transmission | 76° |
| `zemax_kba_R` | Zemax article 42661816933139, attachment `test_example1.zmx` | 0.05 | 5 µm | 550 nm | −45°, 0° from the opposite side | reflection | 22.5° |
| `unslanted_T` | F-PSM note §6.1, case T-a | 0.03 | 10 µm | 500 nm | +33.8627°, −33.8627° | transmission | −90° (unslanted) |

Sweeps (defined in `cases.py`): θ at fixed λ, and λ at the Bragg angle θ = angle of **k₁**.

**What is already established (VERIFIED 2026-09-17).** F-PSM reproduces all 81 numbers of the
note's §6. The Python RCWA reproduces the note's seven printed rigorous efficiencies to ≤ 5×10⁻¹¹,
conserves energy to 10⁻¹³, and is converged (L = 10 vs 14: < 6×10⁻¹³) on all four cases.
OpticStudio's **native Kogelnik** equals the Python Kogelnik to 3.9×10⁻⁸. Against Python RCWA, F-PSM's
rms error is 2.4×10⁻⁵ where Kogelnik's is 1.0×10⁻² (`lumerical_R`). A 1,575-grating stress scan found
four F-PSM failure modes. Full account: `README.md`.

## 2. What is in the pack

```
README.md                       the report so far (read §"What was run, and what could not be")
FPSM_Simple_Implementation.pdf  the F-PSM definition, test vectors, working domain
code/cases.py                   THE grating definitions (single source of truth)
code/fpsm.py kogelnik.py rcwa_te.py          the three engines (proven)
code/test_fpsm_vectors.py test_engines.py    their tests (proven)
code/zos_connect.py             standalone ZOS-API connection, install resolved by scan (proven)
code/zemax_kogelnik_nsc.py      NSC bench + native Kogelnik sweep (proven)
code/zemax_rcwa_nsc.py          RCWA-DLL sweep on the same bench            ** UNTESTED at the DLL **
code/lumerical_rcwa.py          Lumerical RCWA sweep via lumapi             ** UNTESTED draft **
code/compare.py                 figures + summary; overlays/scoring of Ansys RCWA tables (proven with synthetic tables)
code/map2d.py stress_scan.py stress_summary.py    supporting analyses (proven; no need to rerun)
results/                        all existing tables (compare_*, zemax_*, summary.csv, stress_*)
figures/                        all existing figures
zemax/kogelnik_native_*.zmx     the four proven benches, written by OpticStudio
zemax/ansys_kba/*.zar           the Zemax article's archives, as downloaded
zemax/ansys_kba/test_example1/  test_example1 restored (its .zmx references the RCWA DLL)
zemax/ansys_kba/dll/vhg_constant_hologram230313.dll    the RCWA DLL (from the same archive)
lumerical/ansys_example/        the Lumerical article's example: vhg_R_6um_model.fsp, VHG_test_R.zar, licence
```

## 3. What is missing

1. `results/zemax_rcwa_<case>.csv` — OpticStudio RCWA DLL, 4 cases. Blocked before by a licence gate:
   the DLL requires a Premium/Enterprise subscription and **terminates OpticStudio** on a
   Professional seat (no error dialog in headless mode, no Windows fault event).
2. `results/lumerical_rcwa_<case>.csv` — Lumerical RCWA solver, 4 cases. Lumerical was not installed.
3. The overlays, the scores against Ansys RCWA, and the verdict (§5, steps 5–6).

Table schema for both (one schema for every Ansys table; `compare.py` depends on it):
`sweep` (`angle`|`wavelength`), `source_tilt_x_deg` (**= −θ**, degrees), `lambda_um`,
`P_entry_side_W` (reflected first order, per 1 W incident), `P_far_side_W` (transmitted first order).

## 4. Environment

Python ≥ 3.10 (pack built on 3.13.15) with `numpy` (2.2.6), `matplotlib` (3.10.6), `pythonnet`
(3.0.5). `lumapi` ships with Lumerical; `lumerical_rcwa.py` finds it by scan or via the env var
`LUMAPI_PATH` (full path to `lumapi.py`). OpticStudio is found by `zos_connect.py` (env `ZEMAX_ROOT`
= the *install* folder holding `ZOSAPI.dll`, else a scan of `C:\Program Files\*Zemax OpticStudio*`).
Never hard-code an install path. Run every script from inside `code/`.

## 5. Procedure

**Step 0 — prove the pack on this machine.** `python test_fpsm_vectors.py` → `0 failure(s)` (81
checks). `python test_engines.py` → `0 failure(s)` (a `SKIP … PSMGRATING_PATH` line is expected).
Then `python zemax_kogelnik_nsc.py lumerical_R` (~5 min) and `python compare.py lumerical_R`: the
column `zemax_native_vs_my_kogelnik_max` must be ≤ 1×10⁻⁶. This proves ZOS-API, the bench and the
licence on *this* install before anything new is tried. Connection is **standalone** (headless).

**Step 1 — install the DLL.** Copy `zemax/ansys_kba/dll/vhg_constant_hologram230313.dll` into
`<Zemax data folder>\DLL\Diffractive\` (the data folder is `app.ZemaxDataDir`, normally
`Documents\Zemax` — *not* the install folder). Confirm the licence edition first
(`app.LicenseStatus`): anything below Premium will kill the server the moment the DLL loads.
The DLL writes a log `vhg_constant_hologram.txt` beside itself — read it when something is odd.

**Step 2 — probe before sweeping.** `python zemax_rcwa_nsc.py --probe lumerical_R`. It loads the
article's own `test_example1.zmx`, prints every Diffraction-tab parameter (name, reflect, transmit),
then on the index-matched bench reports the first order at Bragg for order +1 and −1, the energy sum
over orders −3…3, and the Bragg efficiency for Max Order 3/5/9/13. Gates in §6 (Z1–Z4). What the
script does and why: it copies **all** of Ansys' shipped parameter values from the template, then
overrides *by name* only the grating (X1…Z2, wavelength, thickness, n1, n2, holo index, dn, holo
type, interpolation = 0). Known facts to lean on: `obj.DiffractionData` exposes `Split`, `DLL`,
`StartOrder`, `StopOrder`, `NumberOfParameters`, `Get/SetReflectParameterValue(i)`,
`GetReflectParameterName(i)`, `GetTransmitParamaterName(i)` *(sic — the API's spelling)*; parameter
indices are **0-based** (VERIFIED 2026-09-17 with `Diff2DSample.dll`); the enum is
`ZOSAPI.Editors.NCE.DiffractionSplitType.SplitByDLL`. The DLL's parameter names, in its own order:
Max Order · X1 Y1 Z1 X2 Y2 Z2 · Const. Wave (um) · Holo Thic (um) · n1 · n2 · Holo index · Index
modulation (dn) · Use Coat File · Shrinkage (unused) · Index Shift (unused) · Rotate Grating ·
Interpolation · Test Mode · Only these orders · Stochastic mode · Layer per period · Unused · Holo
Type. `test_example1.zmx` stores, in that order: 5 · 0 1 −1 0 0 −1 · 0.55 · 5 · 1.5 · 1.5 · 1.5 · 0.05
· 0 · 0 · 0 · 0 · 0 · 1 · 16 · 0 · 10 · 0 · 2 (the name↔value pairing is ASSUMED — the probe's
printout is the authority). The article states the DLL supports **collimated construction beams
only** and recommends Holo Type 2; X1…Z2 are direction-defining points (the example uses unit-scale
numbers, so the script passes P/|P|).

**Step 3 — the four OpticStudio RCWA sweeps.** `python zemax_rcwa_nsc.py <case>` for `lumerical_R`,
`zemax_kba_R`, `lumerical_T`, `unslanted_T` → `results/zemax_rcwa_<case>.csv` and
`zemax/rcwa_dll_<case>.zmx`. For the transmission cases the first-order label may swap past the
fringe-parallel angle (see §7, "order labelling"); handle it as `compare.py` does, do not hide it.

**Step 4 — Lumerical.** `python lumerical_rcwa.py --probe lumerical_R` first: it prints the real
property names of `::model`, `RCWA` and `rectangle` and the result structure. Replace every line
marked `ASSUMED` in `configure()`/`main()` with what the probe and the Lumerical documentation say
(RCWA solver object, `addrcwa`, grating characterisation result). Then run the four cases →
`results/lumerical_rcwa_<case>.csv`, `lumerical/rcwa_<case>.fsp`. `unslanted_T` needs the bypass
described in the script's docstring (K_z = 0 makes the model's z-period infinite). Article settings
to reproduce first, as a sanity anchor: background index 1.5, S polarisation, recording 0° and 28°
in-medium from opposite sides, λ = 550 nm, *n*₀ = 1.5, Dn = 0.04 (the cosine **amplitude**: the
model's script builds `n0 + Dn*cos(...)`), thickness 6 µm by layer repetition, one z-period
(0.1947 µm) sliced into 10 layers, x-period 0.781 µm; the article reports ≈ 0.80 at θ = 0° and −28°.

**Step 5 — overlay and score.** `python compare.py` (all cases). It picks up the new tables
automatically, adds them to `figures/compare_<case>.png` (markers) with their deviation from the
Python RCWA (dotted, lower panels), and appends to `results/summary.csv` for each Ansys solver:
`python_rcwa_vs_<solver>_rms/_max`, `kogelnik_vs_<solver>_rms`, `fpsm_vs_<solver>_rms`.
Look at every figure you regenerate before calling it done.

**Step 6 — write up.** Create `results/ansys_rcwa_verdict.md`: verdict first, in bold; then one
table per solver (case × sweep: Python-RCWA-vs-Ansys rms and max, Kogelnik-vs-Ansys rms,
F-PSM-vs-Ansys rms); then every deviation from this file's procedure, every ASSUMED item you
settled and how, and everything that remains UNTESTED. Update `README.md`: in "What was run, and what
could not be" replace the two **NOT RUN** rows with what you ran (versions, editions, agreement
figures) and adjust the sentence "the RCWA in every figure is the Python referee". Update
`context.md` (State, Acceptance gate, Next) with real state and today's date. Do not touch the
stress-scan text unless you reran it.

**Step 7 — optional, only after 1–6 pass.** The Zemax article's headline picture is an exit-pupil-
expander image, Kogelnik beside RCWA (`zemax/ansys_kba/EPE with RCWA hologram.zar`). Restore it
(`Tools.OpenRestoreZAR()`, `SetFileName`, `SetOutputFolder`, `RunAndWaitForCompletion`), trace it, and
save the RCWA detector image to `figures/epe_rcwa.png`. **F-PSM cannot be added to that picture**:
the waveguide index differs from the hologram index and F-PSM is index-matched only. Say so; do not
force it.

## 6. Acceptance gates — how you know it is right

Reference values from the converged Python RCWA (L = 14), TE, per 1 W incident; "R+1 / T+1" is the
first order on the entry / far face, "T0" the undiffracted beam:

| case | θ | λ | RCWA first order | T0 | Kogelnik | F-PSM |
|---|---|---|---|---|---|---|
| `lumerical_R` | 0° and −28° (Bragg) | 550 nm | R+1 = **0.804947** | 0.1947 | 0.805403 | 0.805021 |
| `lumerical_R` | −14° | 550 | R+1 = **0.023944** | 0.9760 | 0.014014 | 0.023941 |
| `lumerical_R` | 7.5° | 550 | R+1 = **0.100148** | 0.8995 | 0.084730 | 0.100133 |
| `lumerical_T` | 0° (Bragg) | 550 | T+1 = **0.978518** | 0.0139 | 0.987535 | 0.986288 |
| `lumerical_T` | 5° | 550 | T+1 = **0.013852** | 0.9790 | 0.025587 | 0.014074 |
| `zemax_kba_R` | −45° and 0° (Bragg) | 550 | R+1 = **0.873784** | 0.1256 | 0.874565 | 0.873949 |
| `zemax_kba_R` | −22.5° | 550 | R+1 = **0.012827** | 0.9868 | 0.042365 | 0.012826 |
| `unslanted_T` | 33.8627° (Bragg) | 500 | T+1 = **0.586035** | 0.4139 | 0.585763 | 0.586035 |
| `unslanted_T` | 36° | 500 | T+1 = **0.243709** | 0.7562 | 0.253146 | 0.243707 |

The off-Bragg rows are the discriminating ones: a genuine RCWA must land on the RCWA column, **not**
on Kogelnik (−14°: 0.0239 against 0.0140 — a factor 1.7). If an Ansys "RCWA" curve sits on the
Kogelnik column, it is not doing what you think (wrong DLL, interpolation grid, wrong model).

- **Z1** probe lists the 24 parameter names above; template values as listed. *Else:* stop, print what it is.
- **Z2** at Bragg one of order ±1 carries R ≈ 0.805 (`lumerical_R`) on the entry-side detector. *Else:*
  geometry/parameter mapping is wrong — check X1…Z2 sign and Holo Type, the thickness unit, `n1/n2`.
- **Z3** Σ over orders −3…3 of (R + T) = 1.000 ± 1×10⁻³. *Else:* the DLL renormalises when orders are
  restricted, or "Only these orders"/"Stochastic mode" is interfering — then trace all orders with
  Save Rays and separate them by direction (ZRD) or by a Detector Polar instead.
- **Z4** |η(Max Order 9) − η(Max Order 5)| < 1×10⁻⁴ at Bragg. *Else:* raise Max Order for the sweeps
  and say so. Separately raise "Layer per period" (10 → 20 → 40) until the Bragg value moves by
  < 1×10⁻⁴; z-slicing error falls roughly as 1/N².
- **Z5 / L1** full sweeps: `python_rcwa_vs_<solver>_rms` ≤ 1×10⁻³ on every case (both solvers are
  rigorous on the same index-matched problem; the residual should be discretisation only). Expected
  then: `fpsm_vs_<solver>_rms` ≈ `fpsm_rms`, and `kogelnik_vs_<solver>_rms` ≈ 1–2×10⁻².
- **L2** Lumerical converged in harmonics (number of orders / k-vectors) and in layers per z-period
  to the same 1×10⁻⁴ criterion; and the layer-repetition thickness is really *d*: 6 µm is 30.8
  z-periods, so check whether "thickness" mode truncates to whole periods (compare with the article's
  "full structure" set-up, repetition off, 308 layers). A 0.036 µm thickness error moves η by ~10⁻³.
- **If a gate cannot be met,** do not loosen it silently: report the measured value, your diagnosis,
  and mark the affected table `UNTESTED` or `FAILED`.

## 7. Pitfalls already paid for (all VERIFIED 2026-09-17 on OpticStudio 2026 R1 unless marked)

- **Thickness units differ by object.** NSC Hologram Lens `Hologram Thickness` (native Kogelnik) is
  in **lens units (mm)**. The RCWA DLL's parameter is labelled `Holo Thic (um)` — micrometres,
  ASSUMED from the label and from the article's value 5; settle it at gate Z2.
- **NSC minimum relative ray intensity** defaults to 10⁻³ and silently zeroes every efficiency below
  10⁻³. The bench sets `SystemData.NonSequentialData.MinimumRelativeRayIntensity = 1e-12`.
- **Order labelling.** OpticStudio's native model traces, as "order 1", the order with
  **K·k_playback > 0**; on a transmission grating the label swaps to the conjugate order past the
  fringe-parallel angle (θ > 14° for `lumerical_T`). `compare.py` applies the same rule to the Python
  models. The RCWA DLL uses integer orders and may **not** swap — check which integer carries each
  Bragg peak of `lumerical_T` (θ = 0° and θ = 28°) and stitch if necessary, documenting it.
- **Cut-off.** Past cut-off of the first order OpticStudio returns 0; F-PSM returns NaN ("not
  applicable"); the bare Kogelnik formula returns a meaningless number, so `kogelnik.py` returns 0.
- **NCE cell types are per column** (`Diff Order` is Double, `Holo type` Integer): ask
  `cell.DataType`. **pythonnet + ZOS-API hangs at interpreter exit**: end with `os._exit`. Print with
  `flush=True` and redirect to a file, or a killed run loses its output. `DiffractionData.
  GetAvailableDLLs()` hung a headless session once (UNTESTED why) — do not call it.
- `*** FRU__delta_init(): Attempt to start when running!` on stdout is harmless.
- **Bench design** (keep it): source, hologram and both absorbing detectors sit *inside* one
  Rectangular Volume of model index *n*₀, so no ray ever meets a Fresnel interface; Source Ray,
  1 analysis ray, 1 W, `Jx = 1, Jy = 0` (TE for a y–z plane of incidence), polarisation ON.
  Native model: Split NSC Rays OFF, `Diff Order = 1`. RCWA DLL: Split NSC Rays **ON**,
  `StartOrder = StopOrder`.
- **RCWA truncation is geometry-dependent.** L = 10 is ample for these four gratings but is 3.8×10⁻²
  wrong for long-period gratings (slant ≈ 55°); never assume convergence, show it.
- **The Ansys knowledge-base pages return HTTP 403 to plain fetches.** They are readable through
  `https://optics.ansys.com/api/v2/help_center/en-us/articles/<id>.json` with a browser user-agent.
  Article ids: 42661816933139 (Zemax RCWA VHG DLL), 37118679173267 (Lumerical VHG), 42661738883219
  (Kogelnik method in OpticStudio), 42661761082899 (dynamic Lumerical-RCWA ↔ OpticStudio link;
  explains Max Order, Interpolation, Only these orders, Stochastic mode), 42661666095891 (SRG RCWA
  DLLs), 36965778575379 (Lumerical layer repetition).

## 8. Deliverables checklist

- [ ] Step 0 outputs: both test scripts `0 failure(s)`; native-Kogelnik agreement ≤ 1×10⁻⁶ on this PC
- [ ] `results/zemax_rcwa_{lumerical_R,zemax_kba_R,lumerical_T,unslanted_T}.csv` + `zemax/rcwa_dll_*.zmx`
- [ ] `results/lumerical_rcwa_{…same four…}.csv` + `lumerical/rcwa_*.fsp`
- [ ] the probe printouts saved as `results/zemax_rcwa_probe.txt`, `results/lumerical_probe.txt`
- [ ] four regenerated `figures/compare_*.png` showing both Ansys RCWA marker sets; `results/summary.csv` with the eight new columns filled
- [ ] `results/ansys_rcwa_verdict.md`; `README.md` and `context.md` updated
- [ ] the final versions of `zemax_rcwa_nsc.py` and `lumerical_rcwa.py`, with every `ASSUMED`/`UNTESTED` label either removed (and the evidence noted) or still standing and explained
- [ ] (optional) `figures/epe_rcwa.png`

Finish with a short report: verdict; the agreement numbers; what you changed in the two drivers and
why; what is still untested. If something is blocked by a licence, an access control or a missing
product, say exactly what and stop there — do not work around it.
