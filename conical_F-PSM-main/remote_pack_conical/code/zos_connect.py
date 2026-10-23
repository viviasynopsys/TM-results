"""Standalone ZOS-API connection for OpticStudio (pythonnet).

The install root is resolved the way ``0_System/env.ps1`` does it -- the
ZEMAX_ROOT environment variable if set, else the highest-named
``C:\\Program Files\\*Zemax OpticStudio*`` folder that actually holds
ZOSAPI.dll.  Never a literal.  Fails loudly if nothing resolves or the licence
is not valid for the API.
"""
import glob
import os

import clr  # pythonnet


def resolve_zemax_root():
    root = os.environ.get("ZEMAX_ROOT")
    if root and os.path.isfile(os.path.join(root, "ZOSAPI.dll")):
        return root
    cands = [p for p in glob.glob(r"C:\Program Files\*Zemax OpticStudio*")
             if os.path.isfile(os.path.join(p, "ZOSAPI.dll"))]
    if not cands:
        raise RuntimeError("ZEMAX_ROOT did not resolve: no install with ZOSAPI.dll found")
    return sorted(cands)[-1]


class Zos:
    def __init__(self):
        root = resolve_zemax_root()
        self.root = root
        clr.AddReference(os.path.join(root, "ZOSAPI_NetHelper.dll"))
        import ZOSAPI_NetHelper
        if not ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize(root):
            raise RuntimeError(f"ZOSAPI_Initializer.Initialize failed for {root}")
        clr.AddReference(os.path.join(root, "ZOSAPI.dll"))
        clr.AddReference(os.path.join(root, "ZOSAPI_Interfaces.dll"))
        import ZOSAPI
        self.ZOSAPI = ZOSAPI
        self.conn = ZOSAPI.ZOSAPI_Connection()
        self.app = self.conn.CreateNewApplication()
        if self.app is None:
            raise RuntimeError("CreateNewApplication returned None")
        if not self.app.IsValidLicenseForAPI:
            raise RuntimeError(f"Licence not valid for ZOS-API: {self.app.LicenseStatus}")
        self.sys = self.app.PrimarySystem

    def close(self):
        if self.app is not None:
            self.app.CloseApplication()
            self.app = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
