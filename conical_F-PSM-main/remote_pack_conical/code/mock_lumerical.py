"""mock_lumerical -- a stand-in for lumapi.FDTD so the conical driver and the comparison can be
exercised end to end on a PC without Lumerical (dryrun_conical.py).  It is a test double, not a
solver: every result is the conical referee's (rcwa_conical.solve) plus a 1e-6 perturbation, laid
out in the dataset shapes VERIFIED in the 2026-09-18 probe dump, with Lumerical's (theta, phi)
read as the report's (theta, phi) (sign +1) and the x order n equal to the report's order m.
The complex amplitudes follow the ASSUMED convention "first letter = incident polarisation"
(Tsp = s in, p out), flux-normalised; compare_conical.py must detect this, not assume it.
"""
import numpy as np

from rcwa_conical import solve

deg = np.pi / 180.0


class MockFDTD:
    def __init__(self, case, noise=1e-6, seed=0):
        self.case = case
        self.noise = noise
        self.rng = np.random.default_rng(seed)
        dz = 1.947e-7
        self.props = {
            "::model": {"angle 1": 0.0, "angle 2": 28.0, "holo type": "R", "n0": 1.5, "Dn": 0.04, "w_rec": 5.5e-7,
                        "setup script": "shipped script"},
            "rectangle": {"z span": dz, "index": "1.5+0.04*cos(...)"},
            "RCWA": {"z span": 1.5 * dz, "x span": 5e-7, "index": 1.5, "layer repetition": "thickness",
                     "layer repetition start interface": 1, "layer repetition end interface": 2,
                     "layer repetition thickness": np.array([[6e-6]]), "layer repetition multiple": 30.8,
                     "interface reference positions": [["::model::rectangle", "max", 10.0], ["::model::rectangle", "min", 1.0]],
                     "k vectors domain": "rectangular", "max number ku": 3, "max number kv": 0, "max number k vectors": 50,
                     "incident angle": "range", "incident angle table": np.zeros((1, 2)), "minimum theta": -50.0,
                     "maximum theta": 30.0, "theta points": 81, "theta spacing": "linear", "minimum phi": 0.0,
                     "maximum phi": 0.0, "phi points": 1, "wavelength center": 5.5e-7, "wavelength span": 0.0,
                     "use wavelength spacing": 1, "frequency points": 1, "propagation axis": "z",
                     "propagation direction": "both", "mesh cells z": 100},
        }
        self.results = {}
        self.calls = []

    def version(self):
        return "MOCK (no Lumerical here)"

    def switchtolayout(self):
        pass

    def runsetup(self):
        pass

    def save(self, path):
        pass

    def close(self):
        pass

    def getnamed(self, obj, prop=None):
        if prop is None:
            return "\n".join(self.props[obj])
        R = self.props["RCWA"]
        if obj == "RCWA" and prop == "simulation z span":
            rep = float(np.asarray(R["layer repetition thickness"]).ravel()[0])
            return rep + (R["z span"] - self.props["rectangle"]["z span"])
        if obj == "RCWA" and prop == "interface actual positions":
            n = int(R["interface reference positions"][0][2])
            return np.linspace(0.0, 1.0, n + 1)
        if obj == "RCWA" and prop == "incident angles actual":
            if R["incident angle"] == "table":
                return np.asarray(R["incident angle table"], float)
            th = np.linspace(R["minimum theta"], R["maximum theta"], int(R["theta points"]))
            return np.column_stack([th, np.zeros_like(th)])
        return self.props[obj][prop]

    def setnamed(self, obj, prop, val):
        self.props[obj][prop] = val
        self.calls.append((obj, prop))

    def run(self, *args):
        c = self.case
        table = np.asarray(self.getnamed("RCWA", "incident angles actual"), float)
        N = table.shape[0]
        ms = np.arange(-3, 4)
        P = {k: np.zeros((1, N, ms.size, 1, 2)) for k in ("Rs_power", "Ts_power", "Rp_power", "Tp_power")}
        A = {k: np.zeros((1, N, ms.size, 1, 2), dtype=complex) for k in ("Rss", "Rsp", "Rps", "Rpp", "Tss", "Tsp", "Tps", "Tpp")}
        E = {k: np.zeros((1, N, 2)) for k in ("Rs", "Ts", "Rp", "Tp")}
        for j, (th, ph) in enumerate(table):
            R = solve(c["k0"], c["n0"], c["n1"], c["K"], c["Psi"], th * deg, ph * deg, c["d"], M=10)
            g0 = R["gam"][10].real
            for k, m in enumerate(ms):
                i = 10 + m
                nz = self.noise * self.rng.standard_normal(4)
                P["Rs_power"][0, j, k, 0, 0] = R["DE_R"][i, 0] + nz[0]
                P["Ts_power"][0, j, k, 0, 0] = R["DE_T"][i, 0] + nz[1]
                P["Rp_power"][0, j, k, 0, 0] = R["DE_R"][i, 1] + nz[2]
                P["Tp_power"][0, j, k, 0, 0] = R["DE_T"][i, 1] + nz[3]
                fl = np.sqrt(R["gam"][i] / g0) if R["open"][i] else 0.0
                for pre, arr in (("R", R["r"]), ("T", R["t"])):
                    Jm = fl * arr[i]                       # [out, in]
                    A[pre + "ss"][0, j, k, 0, 0] = Jm[0, 0]
                    A[pre + "sp"][0, j, k, 0, 0] = Jm[1, 0]      # s in, p out  (ASSUMED convention)
                    A[pre + "ps"][0, j, k, 0, 0] = Jm[0, 1]
                    A[pre + "pp"][0, j, k, 0, 0] = Jm[1, 1]
            E["Rs"][0, j, 0] = R["DE_R"][:, 0].sum()
            E["Ts"][0, j, 0] = R["DE_T"][:, 0].sum()
            E["Rp"][0, j, 0] = R["DE_R"][:, 1].sum()
            E["Tp"][0, j, 0] = R["DE_T"][:, 1].sum()
        axes = dict(lambda_=np.array([[c["lam_um"] * 1e-6]]), theta=table[:, :1], phi=table[:, 1:], n=ms[:, None].astype(float),
                    m=np.array([[0.0]]), forward=np.array([[1.0], [0.0]]), backward=np.array([[0.0], [1.0]]))
        axes["lambda"] = axes.pop("lambda_")
        self.results["grating_power"] = {**axes, **P}
        self.results["grating_characterization"] = {**axes, **A}
        self.results["total_energy"] = {"lambda": axes["lambda"], "theta": axes["theta"], "phi": axes["phi"],
                                        "forward": axes["forward"], "backward": axes["backward"], **E}

    def getresult(self, obj, name=None):
        if name is None:
            return "\n".join(self.results)
        return self.results[name]
