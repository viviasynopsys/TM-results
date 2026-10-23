"""Second RCWA-DLL probe (written 2026-09-18 after the first probe, results/zemax_rcwa_probe.txt):
 (a) why order 0 reads zero on both detectors (gate Z3): Test Mode and "Only these orders" toggled,
     and all orders -5..5 traced at once;
 (b) TE vs TM (Jx=1 / Jy=1) -- is the DLL polarisation-aware?
 (c) Layer per period 160 / 320 at the three reference angles;
 (d) Holo Thic 5.8413 / 6.0360 um (30 / 31 z-periods) -- the DLL's thickness sensitivity against the
     Python referee's, to test whether the slab is realised as whole periods.
Usage: python zemax_rcwa_probe2.py lumerical_R   (prints; writes nothing but the printout)"""
import os
import sys

from cases import CASES
from zemax_kogelnik_nsc import build, trace
from zemax_rcwa_nsc import configure, diffraction_params, find, read_template
from zos_connect import Zos


def setp(holo, name, v):
    dd = holo.DiffractionData
    live = diffraction_params(dd)
    i = live[find(live, name)][0]
    dd.SetReflectParameterValue(i, float(v))
    dd.SetTransmitParameterValue(i, float(v))


def main(name):
    case = CASES[name]
    z = Zos()
    try:
        template = read_template(z)
        nce, holo, src = build(z, case, volume_hologram=False)
        dd = holo.DiffractionData
        taus = [case["zemax_tilt_bragg_deg"], 14.0, -7.5]
        def at(taus_):
            out = []
            for tau in taus_:
                src.TiltAboutX = float(tau)
                out.append([round(v, 7) for v in trace(z, nce, split=True)])
            src.TiltAboutX = float(case["zemax_tilt_bragg_deg"])
            return out
        print("(a) order 0 / energy", flush=True)
        for tm in (1, 0):
            for oto in (16, 0):
                configure(z, holo, case, template, 0)
                setp(holo, "Test Mode", tm); setp(holo, "Only these orders", oto)
                print(f"   order 0 only, Test Mode={tm} Only these orders={oto}: (R,T)={at(taus)}", flush=True)
                configure(z, holo, case, template, 1)
                setp(holo, "Test Mode", tm); setp(holo, "Only these orders", oto)
                print(f"   order 1 only, Test Mode={tm} Only these orders={oto}: (R,T)={at(taus)}", flush=True)
                configure(z, holo, case, template, 0)
                setp(holo, "Test Mode", tm); setp(holo, "Only these orders", oto)
                dd.StartOrder = -5; dd.StopOrder = 5
                print(f"   orders -5..5,  Test Mode={tm} Only these orders={oto}: (R,T)={at(taus)}", flush=True)
        print("(b) polarisation, Layer per period 10", flush=True)
        configure(z, holo, case, template, 1)
        sd = src.SourcesData
        sd.Jx, sd.Jy = 1.0, 0.0
        print(f"   TE Jx=1: {at(taus)}", flush=True)
        sd.Jx, sd.Jy = 0.0, 1.0
        print(f"   TM Jy=1: {at(taus)}", flush=True)
        sd.RandomPolarization = True
        print(f"   random pol: {at(taus)}", flush=True)
        sd.RandomPolarization = False
        sd.Jx, sd.Jy = 1.0, 0.0
        print("(c) Layer per period 160 / 320", flush=True)
        for lp in (160, 320):
            configure(z, holo, case, template, 1, layers=lp)
            print(f"   Layer per period {lp}: {at(taus)}", flush=True)
        print("(d) thickness sensitivity at Layer per period 40", flush=True)
        for d in (5.8413, 6.0, 6.0360, 6.2):
            configure(z, holo, case, template, 1, layers=40)
            setp(holo, "Thic", d)
            print(f"   Holo Thic {d} um: {at(taus)}", flush=True)
    except Exception:
        import traceback
        traceback.print_exc()
    finally:
        try:
            z.close()
        except Exception as e:
            print("close:", e)
        sys.stdout.flush()
        os._exit(0)


if __name__ == "__main__":
    main(sys.argv[1])
