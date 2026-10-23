"""Step 0 of the remote briefs: is this machine ready?  Prints READY or the first thing missing.
Checks Python >= 3.10, numpy, matplotlib; pythonnet + an OpticStudio install (for the OpticStudio
brief); lumapi (for the Lumerical brief).  Never hard-codes an install path.  2026-10-16."""
import glob
import importlib.util
import os
import sys

ok = True
print("Python", sys.version.split()[0], "(need >= 3.10, Windows, not WSL)")
if sys.version_info < (3, 10):
    ok = False
for mod in ("numpy", "matplotlib"):
    try:
        m = __import__(mod)
        print(f"{mod} {m.__version__}")
    except ImportError:
        print(f"MISSING: {mod}  (pip install {mod})")
        ok = False

which = (sys.argv[1].lower() if len(sys.argv) > 1 else "both")

if which in ("opticstudio", "both"):
    try:
        import clr  # noqa: F401
        print("pythonnet OK")
    except ImportError:
        print("pythonnet MISSING (pip install pythonnet)  -- needed for the OpticStudio brief only")
        ok = ok if which == "both" else False
    root = os.environ.get("ZEMAX_ROOT", "")
    cands = [root] if root else []
    cands += sorted(glob.glob(r"C:\Program Files\*Zemax OpticStudio*"), reverse=True)
    hit = [p for p in cands if p and os.path.isfile(os.path.join(p, "ZOSAPI.dll"))]
    print("OpticStudio install with ZOSAPI.dll:", hit[0] if hit else "NOT FOUND (set ZEMAX_ROOT)")
    if not hit and which == "opticstudio":
        ok = False

if which in ("lumerical", "both"):
    cands = [os.environ.get("LUMAPI_PATH", "")]
    for pat in (r"C:\Program Files\Lumerical\v*\api\python\lumapi.py",
                r"C:\Program Files\ANSYS Inc\v*\Lumerical\api\python\lumapi.py",
                r"C:\Program Files\Ansys Inc\v*\Optical Products\Lumerical\api\python\lumapi.py",
                r"C:\Program Files\ANSYS Inc\v*\Optics\Lumerical\api\python\lumapi.py"):
        cands += sorted(glob.glob(pat), reverse=True)
    hit = [c for c in cands if c and os.path.isfile(c)]
    print("lumapi.py:", hit[0] if hit else "NOT FOUND (set LUMAPI_PATH to its full path)")
    if not hit and which == "lumerical":
        ok = False

here = os.path.dirname(os.path.abspath(__file__))
for f in ("cases.py", "rcwa.py", "kogelnik.py", "fpsm.py", "compare.py"):
    if not os.path.isfile(os.path.join(here, f)):
        print("MISSING:", f)
        ok = False
print("READY" if ok else "NOT READY")
sys.exit(0 if ok else 1)
