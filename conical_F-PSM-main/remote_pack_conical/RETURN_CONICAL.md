# RETURN_CONICAL.md — what to bring back, and how it is merged

*2026-10-23. The pack is a slice of `C:\Users\vivia\Code2\1_Optics\Research\canonical FPSM\` on the
benchmark PC. Everything you produce is merged there; nothing here is rewritten by hand.*

## Bring back (zip the whole pack folder, or at least these)

- `results\lumerical_conical_<case>.csv` for the seven cases (a failed one as `*.FAILED.csv`)
- `results\lumerical_probe.txt` (complete, un-elided)
- the convergence tabulation of step 3, saved as `results\lumerical_conical_converge.txt`
- `lumerical\rcwa_conical_<case>.fsp` as saved by Lumerical
- `results\summary_conical.csv` and `figures\conical_<case>.png` (regenerated on the benchmark PC
  from your tables; yours are checked, not kept)
- `results\lumerical_conical_verdict.md`
- the final `code\lumerical_rcwa_conical.py` and, if touched, `code\lumerical_rcwa.py` (every
  `ASSUMED`/`UNTESTED` label either confirmed with the evidence noted, or still standing and explained)
- a short final report in the chat: verdict; the C2 and C5 numbers; what you changed and why; what
  is still untested; Lumerical version and licence edition

## How it is merged on the benchmark PC

1. Your tables, probe, convergence printout and verdict are copied into `canonical FPSM\results\`
   (Ansys tables are never edited; a table that fails its gate stays as `*.FAILED.csv` beside the
   verdict that explains it).
2. `python compare_conical.py` is run there; the scores land in `results\summary_conical.csv` and
   the figures are regenerated.
3. Your driver version is diffed against the pack's; every confirmed label is carried over with
   your date.
4. `REPORT_CONICAL.md` is extended with the Ansys leg and the findings go to the model's author.

Do not rename files, do not change the CSV schema, do not reorder or drop cases.
