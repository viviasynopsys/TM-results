# CLAUDE.md — read this first

You are Claude (Sonnet 5.5, Claude Code) on a PC that has Ansys Lumerical installed. This folder
is a self-contained task pack from a larger research project. You bring results back; you do not
rewrite the project.

1. Read `README_PACK.md` (two minutes), then the whole of `EXECUTION_CONICAL_LUMERICAL.md` before
   running anything. Run every script from inside `code/`, with a Windows Python (not WSL).
2. Rules that apply to everything you write: British English; every claim carries a status and a
   date (`VERIFIED <date>`, `ASSUMED`, `UNTESTED`); never tune or fabricate a number; fail loudly
   (a failed run writes `*.FAILED.csv`, never overwrites a good table); never hand-edit `.fsp` or
   any binary; verdict first, evidence after; do not modify anything under `reference/`,
   `lumerical/ansys_example/`, `results/reference_conical.csv` or `evidence/`.
3. If a licence, an access control or a missing product blocks a step, say exactly what and stop
   there. Do not work around a licence.
4. When done, follow `RETURN_CONICAL.md`.
