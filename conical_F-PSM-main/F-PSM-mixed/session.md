# session.md — mixed phase–absorption leg on the Lumerical PC

*Session 2026-10-23 23:35 to 2026-10-24 02:15 (BST), written 2026-10-24. Solver PC: Windows, Ansys Lumerical 2026 R1 (build 8.35.4609), Python 3.13.14 (numpy 2.5.1, scipy 1.18.1, matplotlib 3.11.1). British English. Every claim carries a status: VERIFIED `<date>`, ASSUMED or UNTESTED.*

## 1. The request and the outcome

Request: run the pack in `remote_pack_mixed`, reading `RESTORE_EXTENSIONS.md` first and following it. Chain followed: `RESTORE_EXTENSIONS.md` → `CLAUDE.md` → `README_PACK.md` → `EXECUTION_MIXED_LUMERICAL.md` → `RETURN_MIXED.md`.

Outcome (VERIFIED 2026-10-24): the brief was carried out through step 4 for sets A, B and C (58 configs, no failures) and the verdict is written (`results/lumerical_mixed_verdict.md`). **Lumerical's RCWA solves an absorbing index through an `import` primitive and agrees with the independent complex-permittivity conical referee: rms ≤ 5.1 × 10⁻⁵ on the 54 grating configs of sets A and B, ≤ 2.1 × 10⁻⁴ on set C; gates X0–X4 and X6 pass.** The `equation` route is rejected. Set D (optional) and gate X5 were not run: UNTESTED.

## 2. What was done, in order

| time | step | command or action | result |
|---|---|---|---|
| 23:35 | restore | `restore.ps1` (saved from `RESTORE_EXTENSIONS.md`; no Python on PATH for the first attempt, `py` found 3.13 afterwards) | `restored 94 OK, 0 BAD` |
| 23:37–23:59 | step 0 | `check_env.py lumerical` (with `LUMAPI_PATH` set), `test_mixed.py` and `test_conical.py` (run concurrently, so slower than the brief's 1 and 6 minutes), `lumerical_rcwa_mixed.py --list`, `compare_mixed.py` | READY; 0 failures of 11 gates (665 s); 0 failures of 12 gates (1287 s); 100 configs; "no Lumerical table yet" |
| 23:59 | step 1a | `--probe-loss` (route `equation`) | **FAILED**: all three tokens × two signs; the equation reads back real, the imaginary part is not stored |
| 00:00 | step 1b | `--probe-loss --index-route import` | P01 passed (1.44 × 10⁻¹³ against the Fresnel oracle), then crashed in `remove_import` (analysis mode). Driver fixed (section 4), re-run: P02 passed (6.72 × 10⁻⁵). Wrote `loss_route.json` |
| 00:04 | step 2 | `--converge A05_…` | X6 passes (section 3) |
| 00:05–02:10 | step 3 | `--set A`, `--set B`, `--set C` in one chain | 30 + 24 + 4 tables, 0 `*.FAILED.csv`; `.fsp` per set saved by Lumerical |
| 02:11 | step 4 | `compare_mixed.py` | `summary_mixed.csv`, 58 figures; gates scored with a pandas script (kept in `session_tools/` of the copy) |
| 02:15 | step 5 | verdict written | `results/lumerical_mixed_verdict.md` |

Environment used for every solver command: `LUMAPI_PATH = C:\Program Files\ANSYS Inc\v261\Lumerical\api\python\lumapi.py`, run from `code\`.

## 3. Results in brief

| item | value |
|---|---|
| route / sign / token | `import` / +1 (`n + iκ` absorbs) / `n+ik array` (VERIFIED 2026-10-23) |
| P1 | max \|T0 − Fresnel oracle\| = 1.44 × 10⁻¹³ (gate 10⁻³) |
| P2 | R+1 TE 0.735904, TM 0.668684; power sums 0.915632 / 0.909348; max deviation 6.72 × 10⁻⁵ (gate 2 × 10⁻³) |
| X0 | `power_max` ≤ 0.9293; absorbed-fraction rms ≤ 4.6 × 10⁻⁵ |
| X1 | rms: A ≤ 4.3 × 10⁻⁵ (s) / 5.1 × 10⁻⁵ (p); B ≤ 7 × 10⁻⁶ / 9 × 10⁻⁶; C ≤ 1.84 × 10⁻⁴ / 2.10 × 10⁻⁴ |
| X2 | amplitude convention "in-first, power-normalised" on 58; time convention `exp(−iωt)` conjugated on 57 (the other one is a zero-signal config); Jones error ≤ 3.55 × 10⁻⁴ |
| X3 | `jones_within_bound` True on 58 of 58 |
| X4 | `within_E19` True on 30 of 30; observed/E19 in [0.99984, 1.0000014] |
| X5 | UNTESTED (set D not run) |
| X6 | harmonics 3 → 5: no change; layers 80/160/320: η_s steps 1.47 × 10⁻⁵ then 3.5 × 10⁻⁶ (ratio 4.2) |
| δ orientation | slanted T −1, slanted R +1, unslanted T −1 (`results/delta_orientation.json`) |
| one-sided configs (A10, A20, A30) | signal ≤ 5.5 × 10⁻¹³ (brief: ≤ 10⁻⁶) |
| observation | unslanted gratings agree to 10⁻⁷–10⁻⁵, not the brief's 10⁻¹¹; cause UNTESTED (likely the 256-point import grid) |

## 4. Changes to the pack

1. `code/lumerical_rcwa_mixed.py` (the only file of the shipped 94 that differs from its SHA-256 in the restore table, by design):
   - `remove_import()` now calls `fdtd.switchtolayout()` before selecting and deleting the import object; the session is in analysis mode after a solve (VERIFIED 2026-10-23: the first import-route probe stopped on `you cannot modify simulation objects while in analysis mode`).
   - docstring and inline labels updated: the `UNTESTED`/`ASSUMED` on the complex index, the equation tokens and `importnk2` replaced by the evidence and date. No logic change otherwise.
2. `restore.ps1` added at the pack root (the script from `RESTORE_EXTENSIONS.md`, unchanged). A re-run now reports `93 OK, 1 BAD` (the edited driver) and nothing else.
3. New: `results/` (tables, probe logs, route and δ files, convergence printout, summary, verdict), `figures/` (58 PNG), `lumerical/rcwa_mixed_{A,B,C}.fsp` (saved by Lumerical, not hand-edited), this file.

## 5. Incident: `evidence/` was overwritten, then restored

`CLAUDE.md` forbids modifying `evidence/`. Step 0 of the brief runs `test_conical.py` and `test_mixed.py`, and both write `../evidence/conical_gates.txt` and `../evidence/mixed_gates.txt` themselves (stated in their headers). They replaced the two shipped files with this PC's output (different SHA-256: timings differ). It was noticed on 2026-10-24 when listing the files changed in the session.

Repair (VERIFIED 2026-10-24): the two original files were re-extracted from `Downloads\remote_pack_mixed_2026-10-23_txt.zip`, their SHA-256 checked against the restore table (`7afbee97…` and `9ec87299…`), and written back to `evidence/`. This PC's own gate output, which is genuine evidence of the run (0 failures of 11 and of 12), is kept as `evidence_solver_pc/conical_gates_solver_pc_2026-10-23.txt` and `mixed_gates_solver_pc_2026-10-23.txt` in `F-PSM-mixed`. For any future run on a pack copy: copy `evidence/` aside first, or run those two scripts with the output redirected.

## 6. Files this session produced

All under `remote_pack_mixed\` unless noted. A copy of everything below, except the unchanged restore chain, is in `..\F-PSM-mixed\` (same layout).

- `session.md` (this file), `restore.ps1`
- `code\lumerical_rcwa_mixed.py` (modified)
- `results\lumerical_mixed_<case_id>.csv` × 58 (A01–A30, B01–B24, C01–C04), `summary_mixed.csv`, `loss_route.json`, `delta_orientation.json`, `lumerical_mixed_probe_loss.txt`, `lumerical_mixed_probe_loss_equation.FAILED.txt`, `lumerical_mixed_converge_A05_slanted_T_n024_k01200_01500_d000_phi150.txt`, `lumerical_mixed_verdict.md`
- `figures\mixed_<case_id>.png` × 58 (no `mixed_phase_scan_*` figures: set D not run)
- `lumerical\rcwa_mixed_A.fsp`, `rcwa_mixed_B.fsp`, `rcwa_mixed_C.fsp`
- only in the copy: `evidence_solver_pc\` (this PC's two gate outputs) and `session_tools\` (`gates.py`, the X0–X4 checker; `sheets.py`, the contact-sheet maker; `sheet_1..5.png`, the figure contact sheets used for the figure review)

Not kept: the console logs of the runs (temporary files, deleted); every number in them is in the CSVs, `loss_route.json`, the convergence printout and the verdict. The traceback of the first import-route probe is quoted in section 2 only.

## 7. What remains UNTESTED

- Set D (40 configs): D01 θ sweep, the X5 phase scan, and the δ orientation of unslanted T at φ = 135°.
- Convergence of the import grid (256 × 2 × 1024).
- The Lumerical licence edition (not read out; the RCWA solver ran with no licence dialog or error).
- κ above 0.02 and thickness above 13.3 µm.

## 8. To resume or extend

```
cd code
set LUMAPI_PATH=C:\Program Files\ANSYS Inc\v261\Lumerical\api\python\lumapi.py
py lumerical_rcwa_mixed.py --set D          # optional, about 25–30 min; good tables are skipped, --redo forces
py compare_mixed.py                         # rescores everything, adds mixed_phase_scan_<family>.png
```

`results\loss_route.json` is read by every run, so no route flags are needed. To return the work, follow `RETURN_MIXED.md` (zip the whole pack folder).
