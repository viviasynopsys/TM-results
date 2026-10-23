# README_PACK — the conical (out-of-plane) leg of the F-PSM volume-grating benchmark

*Pack built 2026-10-23 on the benchmark PC from `1_Optics/Research/canonical FPSM/`. British English.*

## What this is about, in one paragraph

A closed-form diffraction-efficiency model for volume holographic gratings, **F-PSM**, was
benchmarked in September and October 2026 for in-plane replay (TE, then TM) against Kogelnik's
theory, an independent Python RCWA referee and the RCWA solvers of Ansys Lumerical and Ansys
Zemax OpticStudio. On 20 October 2026 the model's author extended F-PSM to the **full conical
geometry** (replay out of the plane of the grating vector: polar angle θ and azimuth φ, both
polarisations coupled) and sent two technical reports with the equations. On the benchmark PC
those equations were transcribed, an independent conical vector RCWA referee was written, and the
author's own 504-point sweep and every certified bound in his reports were reproduced to four
figures (`evidence/conical_gates.txt`, 12 gates, VERIFIED 2026-10-23). **What is missing is the
Ansys data**: Lumerical's RCWA solver at conical incidence on the same gratings. This pack exists
to produce it. Nothing in this pack asks you to invent a formula or tune a number.

## The brief

| product on this machine | brief | produces |
|---|---|---|
| Ansys Lumerical (2025 R1 or later, RCWA solver licensed) | `EXECUTION_CONICAL_LUMERICAL.md` | `results/lumerical_conical_<case>.csv` ×7, a probe dump, a convergence check, `results/summary_conical.csv`, the figures, a verdict |

## Layout

```
CLAUDE.md                          start here
README_PACK.md                     this file
EXECUTION_CONICAL_LUMERICAL.md     the brief: procedure, gates, the numbers to hit, pitfalls
RETURN_CONICAL.md                  what to bring back and how it is merged
code/                              every script; run from here
  cases_conical.py                 THE case definitions (single source of truth): four benchmark
                                   gratings on their conical Bragg cones + the author's three families
  rcwa_conical.py  basis.py        the conical vector RCWA referee (gated, VERIFIED 2026-10-23)
  fpsm_conical.py                  the conical F-PSM closed forms, leading and corrected (transcribed, gated)
  test_conical.py  sweep_david.py  the gates (run test_conical.py in step 0)
  lumerical_rcwa_conical.py        the conical Lumerical driver      ** UNTESTED at the solver **
  lumerical_rcwa.py                the September/October planar driver it builds on (VERIFIED 2026-09-18)
  compare_conical.py               scores + figures; run after every case
  mock_lumerical.py  dryrun_conical.py   the test double used on the benchmark PC; not needed here
  cases.py rcwa.py fpsm.py kogelnik*.py test_*.py ...   the planar benchmark code (unchanged)
results/reference_conical.csv      THE NUMBERS TO HIT: referee + closed forms on every row you will solve
results/sweep_504.csv              the author's sweep re-run here (for reference only)
evidence/conical_gates.txt         the 12 gates passed on the benchmark PC
reference/                         the two reports (PDF), the September verdicts and article texts,
                                   the 2026-09-18 probe dump (every property and result name)
lumerical/ansys_example/           the Lumerical article's model (vhg_R_6um_model.fsp) -- do not modify
figures/                           the benchmark PC's figures (yours are added by compare_conical.py)
```
