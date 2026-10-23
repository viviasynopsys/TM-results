"""dryrun_conical -- exercise lumerical_rcwa_conical.py and compare_conical.py end to end with the
mock solver (no Lumerical needed).  Outputs go to ../evidence/dryrun/, never to results/, so they can
never be mistaken for solver tables.

    python dryrun_conical.py
Pass criteria: every case runs to a table, the planar sign fit returns +1 / order +1, the s and p
rms against the referee are at the mock's 1e-6 noise, and the amplitude convention is detected as
"in-first, power".
"""
import os
import shutil
import sys

import lumerical_rcwa as sept
import lumerical_rcwa_conical as drv
import compare_conical as cmp
from mock_lumerical import MockFDTD

EV = os.path.join(drv.ROOT, "evidence", "dryrun")
os.makedirs(EV, exist_ok=True)
drv.RES = EV
cmp.RES = EV
cmp.FIG = EV
shutil.copy(os.path.join(drv.ROOT, "results", "reference_conical.csv"), EV)

cases = drv.all_cases()
names = sys.argv[1:] or ["lumerical_R", "unslanted_T", "Slanted_R"]
sept.load_lumapi = lambda: None
for name in names:
    sept.open_model = lambda lumapi, path, _c=cases[name][0]: MockFDTD(_c)
    print(f"=== dry run {name} ===", flush=True)
    drv.main([name, "--harmonics", "3", "--layers", "160"])
print("=== compare ===", flush=True)
cmp.main(names)
