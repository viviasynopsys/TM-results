"""build_pack -- assemble remote_pack_conical/ for the Lumerical PC from the 2026-10-16 TM pack
(verified driver, model, article texts, probe dump) plus everything conical built here, and zip it.

    python build_pack.py            # (re)build ../remote_pack_conical/ ; docs come from ../pack_docs/
    python build_pack.py --zip      # also write the two archives and verify the second one:
        ../remote_pack_conical_2026-10-23.zip        the pack as is
        ../remote_pack_conical_2026-10-23_txt.zip    every file renamed <name>.txt, plus
                                                     RESTORE_EXTENSIONS.md at the root: the manifest
                                                     (file, original name, original extension,
                                                     SHA-256) and the agent instructions and scripts
                                                     that restore the names.  The builder extracts it
                                                     to a temporary folder, runs the shipped Python
                                                     restore script and checks the result against the
                                                     pack file by file before reporting success.
Stray folders (__pycache__, evidence/dryrun) are never shipped.
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TM = os.path.normpath(os.path.join(ROOT, "..", "2026-10-16", "remote_pack"))
PACK = os.path.join(ROOT, "remote_pack_conical")
DATE = "2026-10-23"
NEW_CODE = ["basis.py", "rcwa_conical.py", "fpsm_conical.py", "sweep_david.py", "test_conical.py", "cases_conical.py",
            "lumerical_rcwa_conical.py", "compare_conical.py", "mock_lumerical.py", "dryrun_conical.py",
            "figures_conical.py", "build_pack.py"]
DOCS = ["CLAUDE.md", "README_PACK.md", "EXECUTION_CONICAL_LUMERICAL.md", "RETURN_CONICAL.md"]
SKIP_DIRS = ("__pycache__", "dryrun")
MANIFEST = "RESTORE_EXTENSIONS.md"
ZIP_PLAIN = os.path.join(ROOT, f"remote_pack_conical_{DATE}.zip")
ZIP_TXT = os.path.join(ROOT, f"remote_pack_conical_{DATE}_txt.zip")


def copytree(src, dst):
    for dp, dn, fn in os.walk(src):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        rel = os.path.relpath(dp, src)
        out = os.path.join(dst, rel) if rel != "." else dst
        os.makedirs(out, exist_ok=True)
        for f in fn:
            shutil.copy2(os.path.join(dp, f), os.path.join(out, f))


def build():
    os.makedirs(PACK, exist_ok=True)
    for sub in ("code", "results", "evidence", "figures", "lumerical", "reference"):
        os.makedirs(os.path.join(PACK, sub), exist_ok=True)
    for f in os.listdir(os.path.join(TM, "code")):
        if f.endswith(".py"):
            shutil.copy2(os.path.join(TM, "code", f), os.path.join(PACK, "code", f))
    copytree(os.path.join(TM, "lumerical"), os.path.join(PACK, "lumerical"))
    copytree(os.path.join(TM, "reference"), os.path.join(PACK, "reference"))
    for f in ("tm_reference_values.md", "tm_referee_gates.txt", "kogelnik_gates.txt"):
        p = os.path.join(TM, "evidence", f)
        if os.path.isfile(p):
            shutil.copy2(p, os.path.join(PACK, "reference", "tm_pack_" + f))
    for f in NEW_CODE:
        shutil.copy2(os.path.join(HERE, f), os.path.join(PACK, "code", f))
    for f in ("reference_conical.csv", "sweep_504.csv", "mg_fig5_thickness_scan.csv"):
        shutil.copy2(os.path.join(ROOT, "results", f), os.path.join(PACK, "results", f))
    shutil.copy2(os.path.join(ROOT, "evidence", "conical_gates.txt"), os.path.join(PACK, "evidence", "conical_gates.txt"))
    for f in ("sweep_errors_vs_mu.png", "reference_conical_bench.png"):
        shutil.copy2(os.path.join(ROOT, "figures", f), os.path.join(PACK, "figures", f))
    for f in ("Conical_F_PSM_RCWA_Report.pdf", "Conical_F_PSM_Mixed_Report.pdf"):
        shutil.copy2(os.path.join(ROOT, f), os.path.join(PACK, "reference", f))
    for d in DOCS:
        shutil.copy2(os.path.join(ROOT, "pack_docs", d), os.path.join(PACK, d))
    # strays left by runs made inside the pack (gates, dry run) are removed, never shipped
    for dp, dn, fn in os.walk(PACK):
        for d in list(dn):
            if d in SKIP_DIRS:
                shutil.rmtree(os.path.join(dp, d))
                dn.remove(d)
    print(f"pack: {len(pack_files())} files in {PACK}")


def pack_files():
    """[(absolute path, pack-relative path with forward slashes)], sorted, strays excluded."""
    out = []
    for dp, dn, fn in os.walk(PACK):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        for f in fn:
            if f.endswith(".pyc"):
                continue
            p = os.path.join(dp, f)
            out.append((p, os.path.relpath(p, PACK).replace(os.sep, "/")))
    return sorted(out, key=lambda t: t[1])


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


ROW_RE = r"^\| `([^`]+)` \| `([^`]+)` \| `[^`]*` \| `([0-9a-f]{64})` \|"

RESTORE_PY = f'''import os, hashlib, re
ok = bad = 0
for line in open("{MANIFEST}", encoding="utf-8"):
    m = re.match(r"{ROW_RE}", line)
    if not m:
        continue
    src, dst, want = m.groups()
    if os.path.isfile(src):
        os.replace(src, dst)
    if not os.path.isfile(dst):
        print("MISSING", dst); bad += 1; continue
    got = hashlib.sha256(open(dst, "rb").read()).hexdigest()
    if got == want:
        ok += 1; print("OK ", dst)
    else:
        bad += 1; print("BAD", dst)
print(f"restored {{ok}} OK, {{bad}} BAD")
'''

RESTORE_PS = f'''$ok = 0; $bad = 0
Get-Content {MANIFEST} -Encoding UTF8 | ForEach-Object {{
  if ($_ -match '{ROW_RE}') {{
    $src = $Matches[1]; $dst = $Matches[2]; $want = $Matches[3]
    if (Test-Path -LiteralPath $src) {{ Move-Item -LiteralPath $src -Destination $dst -Force }}
    if (-not (Test-Path -LiteralPath $dst)) {{ Write-Output "MISSING $dst"; $bad++; return }}
    $got = (Get-FileHash -LiteralPath $dst -Algorithm SHA256).Hash.ToLower()
    if ($got -eq $want) {{ $ok++; Write-Output "OK  $dst" }} else {{ $bad++; Write-Output "BAD $dst" }}
  }}
}}
Write-Output "restored $ok OK, $bad BAD"
'''


def manifest_text(files):
    n = len(files)
    lines = [
        f"# {MANIFEST} - read this first, then restore the file names",
        "",
        f"*Written {DATE} on the benchmark PC. Every file in this archive except this one has had `.txt` appended to",
        "its name so that e-mail accepts the archive. Nothing else was changed: the SHA-256 in the table is of each",
        "file's contents, which renaming does not touch. British English.*",
        "",
        "## Agent instructions",
        "",
        "1. Work inside the extracted `remote_pack_conical` folder; this file sits at its root.",
        "2. Run ONE of the two scripts below from that folder. Each reads the table at the end of this file, renames",
        "   every file back to its original name, checks its SHA-256, prints one line per file and a final count.",
        "   Running it twice is harmless: a file already restored is only re-checked.",
        f"3. The run must end with `restored {n} OK, 0 BAD`. If any file is `BAD` or `MISSING`, stop and report it;",
        "   do not continue with a damaged pack.",
        "4. Then open `CLAUDE.md` (now restored) and follow it.",
        "",
        "### Python (3.8 or later)",
        "",
        "Save as `restore.py` beside this file and run `python restore.py`:",
        "",
        "```python",
        RESTORE_PY.rstrip(),
        "```",
        "",
        "### PowerShell",
        "",
        "Save as `restore.ps1` beside this file and run `powershell -ExecutionPolicy Bypass -File restore.ps1`:",
        "",
        "```powershell",
        RESTORE_PS.rstrip(),
        "```",
        "",
        f"## Table: {n} files",
        "",
        "Paths are relative to the pack root, forward slashes. `original extension` is informational; the",
        "scripts use `original name`.",
        "",
        "| file in this archive | original name | original extension | sha256 |",
        "|---|---|---|---|",
    ]
    for p, rel in files:
        ext = os.path.splitext(rel)[1] or "(none)"
        lines.append(f"| `{rel}.txt` | `{rel}` | `{ext}` | `{sha(p)}` |")
    return "\n".join(lines) + "\n"


def zips():
    files = pack_files()
    with zipfile.ZipFile(ZIP_PLAIN, "w", zipfile.ZIP_DEFLATED) as z:
        for p, rel in files:
            z.write(p, f"remote_pack_conical/{rel}")
    print(f"-> {ZIP_PLAIN}  ({len(files)} files)  sha256 {sha(ZIP_PLAIN)}")
    with zipfile.ZipFile(ZIP_TXT, "w", zipfile.ZIP_DEFLATED) as z:
        for p, rel in files:
            z.write(p, f"remote_pack_conical/{rel}.txt")
        z.writestr(f"remote_pack_conical/{MANIFEST}", manifest_text(files))
    print(f"-> {ZIP_TXT}  ({len(files)} files + {MANIFEST})  sha256 {sha(ZIP_TXT)}")


def verify_txt():
    """Extract the .txt archive into a temporary folder, run the SHIPPED Python restore script there,
    and compare the restored tree with the pack file by file."""
    files = pack_files()
    want = {rel: sha(p) for p, rel in files}
    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(ZIP_TXT) as z:
            z.extractall(tmp)
        root = os.path.join(tmp, "remote_pack_conical")
        with open(os.path.join(root, "restore.py"), "w", encoding="utf-8") as f:
            f.write(RESTORE_PY)
        r = subprocess.run([sys.executable, "-I", "restore.py"], cwd=root, capture_output=True, text=True)
        last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()
        got = {}
        for dp, dn, fn in os.walk(root):
            for f in fn:
                rel = os.path.relpath(os.path.join(dp, f), root).replace(os.sep, "/")
                if rel in ("restore.py", MANIFEST):
                    continue
                got[rel] = sha(os.path.join(dp, f))
        same = got == want
        print(f"verify: restore script says '{last}'; restored tree {'IDENTICAL to' if same else 'DIFFERS from'} the pack "
              f"({len(got)} files restored, {len(want)} expected)")
        if not same:
            for k in sorted(set(got) ^ set(want)):
                print("   only on one side:", k)
            for k in sorted(set(got) & set(want)):
                if got[k] != want[k]:
                    print("   hash differs:", k)
            raise SystemExit("round trip FAILED")
        if last != f"restored {len(want)} OK, 0 BAD":
            raise SystemExit("round trip FAILED: " + last)


if __name__ == "__main__":
    build()
    if "--zip" in sys.argv:
        zips()
        verify_txt()
