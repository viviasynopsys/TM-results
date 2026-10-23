"""Print the stress-scan verdict table from ../results/stress_scan.csv (and save
it as ../results/stress_summary.txt).  Gratings are grouped by their ACTUAL
geometry at Bragg (F-PSM's own eq. 8 rule), not by the psi range they came from."""
import csv, io, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "stress_scan.csv"))))
F = lambda r, k: float(r[k]) if r.get(k, "") not in ("", "nan") else np.nan
out = io.StringIO()
P = lambda *a: print(*a, file=out)

unref = [r for r in rows if r["geometry_at_bragg"] == "unrefereed"]
P(f"gratings: {len(rows)}   unrefereed (RCWA ladder 10..54 not converged to 1e-8): {len(unref)}")
for r in unref:
    P(f"   unrefereed: psi={r['psi_deg']} n1={r['n1']} d={r['d_um']}um residual={r['rcwa_residual']}")
Ls = {}
for r in rows:
    Ls[r["rcwa_L"]] = Ls.get(r["rcwa_L"], 0) + 1
P("RCWA truncation used (L: count):", dict(sorted(Ls.items(), key=lambda kv: int(kv[0]))))

med = lambda L, k: float(np.nanmedian([F(r, k) for r in L])) if L else np.nan
mx = lambda L, k: float(np.nanmax([F(r, k) for r in L])) if L else np.nan
for geo in ("reflection", "transmission"):
    rr = [r for r in rows if r["geometry_at_bragg"] == geo and int(r["n_points"]) > 0]
    P(f"\n=== {geo} at Bragg: {len(rr)} gratings")
    P(f"  gratings with >=1 closed-channel band point: {sum(int(r['n_closed']) > 0 for r in rr)};"
      f" with T/R flag flips inside the band: {sum(int(r['n_flag_flips']) > 0 for r in rr)}")
    P(f"  F-PSM eta outside [0,1] by > 1e-6: {sum(F(r, 'fpsm_out_of_unit') > 1e-6 for r in rr)} gratings, worst excess {mx(rr, 'fpsm_out_of_unit'):.2e}")
    worse = [r for r in rr if F(r, "fpsm_rms") > F(r, "kogelnik_rms")]
    P(f"  F-PSM worse than Kogelnik (whole band): {len(worse)} / {len(rr)}; of those with the Bragg point in-domain: {sum(int(r['bragg_in_domain']) for r in worse)}")
    ind = [r for r in rr if int(r["n_in_domain"]) > 0]
    P(f"  in-domain POINTS only ({sum(int(r['n_in_domain']) for r in rr)} points in {len(ind)} gratings): F-PSM rms median {med(ind, 'fpsm_rms_in_domain_pts'):.1e},"
      f" worst grating {mx(ind, 'fpsm_rms_in_domain_pts'):.1e}, worst point {mx(ind, 'fpsm_max_in_domain_pts'):.1e};  Kogelnik median {med(ind, 'kogelnik_rms_in_domain_pts'):.1e}")
    P("  by modulation (whole band, in and out of domain):")
    P("    n1      F-PSM median   F-PSM max    Kogelnik median   F-PSM worse than Kogelnik")
    for n1 in sorted({float(r["n1"]) for r in rr}):
        L = [r for r in rr if float(r["n1"]) == n1]
        P(f"    {n1:<7g} {med(L, 'fpsm_rms'):.1e}        {mx(L, 'fpsm_rms'):.1e}      {med(L, 'kogelnik_rms'):.1e}           {sum(F(r, 'fpsm_rms') > F(r, 'kogelnik_rms') for r in L)}/{len(L)}")
    P("  by thickness:")
    for d in sorted({int(r["d_um"]) for r in rr}):
        L = [r for r in rr if int(r["d_um"]) == d]
        P(f"    d={d:<3d}um F-PSM median {med(L, 'fpsm_rms'):.1e}  Kogelnik median {med(L, 'kogelnik_rms'):.1e}")
    P("  by slant (median over n1, d):")
    for p in sorted({int(r["psi_deg"]) for r in rr}):
        L = [r for r in rr if int(r["psi_deg"]) == p]
        P(f"    psi={p:<5d} min(gamma)/beta at Bragg {F(L[0], 'bragg_min_gamma_over_beta'):.3f}   F-PSM median {med(L, 'fpsm_rms'):.1e} max {mx(L, 'fpsm_rms'):.1e}   Kogelnik median {med(L, 'kogelnik_rms'):.1e}")
    P("  worst 5 by F-PSM rms:")
    for r in sorted(rr, key=lambda r: -F(r, "fpsm_rms"))[:5]:
        P(f"    psi={r['psi_deg']:>4} n1={r['n1']:<6} d={r['d_um']:>2}um  F-PSM {F(r, 'fpsm_rms'):.2e}  Kogelnik {F(r, 'kogelnik_rms'):.2e}  failmask(T6..T1|R6..R1)={int(r['domain_fail_or']):06b}")
text = out.getvalue()
sys.stdout.write(text)
open(os.path.join(ROOT, "results", "stress_summary.txt"), "w", encoding="utf-8").write(text)
