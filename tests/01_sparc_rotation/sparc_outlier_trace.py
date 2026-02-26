#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - SPARC OUTLIER TRACE AUDIT V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Multiplicative Boost verification against high-mass under-prediction.
# Protocol: Targeted Trace Capture for NGC 2841, 5005, and 3198.
# =============================================================================

import numpy as np
import pandas as pd
import io, os, sys

# 1. SYSTEM LINKING
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
try:
    from network_g_core import engineV5 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: EngineV5 linked for Outlier Audit.")
except ImportError:
    print("❌ ERROR: Physics engine missing.")
    sys.exit(1)

# CONFIG: THE BIG THREE TARGETS
OUTLIERS = ["NGC2841", "NGC5005", "NGC3198"]
SIGMA_STRUCT, SIGMA_INT = 11.5, 3.0
UPS_GRID = np.linspace(0.4, 1.0, 25) 

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")

def load_sparc_mrt(path, columns, skip_lines):
    if not os.path.exists(path): return pd.DataFrame(columns=columns)
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln for ln in f.readlines()[skip_lines:] if ln.strip() and not ln.strip().startswith(('#', '|'))]
    return pd.read_csv(io.StringIO("".join(lines)), sep=r"\s+", names=columns, engine="python", na_values=["----", "...", "NaN"])

if __name__ == "__main__":
    meta = load_sparc_mrt(os.path.join(DATA_DIR, "sparc_metadata.csv"), ['Galaxy','T','D','e_D','f_D','Inc','e_Inc','L36','e_L36','Reff','SBeff','Rdisk','SBdisk','MHI','RHI','Vflat','e_Vflat','Q','Ref'], 34)
    mm = load_sparc_mrt(os.path.join(DATA_DIR, "sparc_mass_models.csv"), ['Galaxy','D_mm','R','Vobs','errV','Vgas','Vdisk','Vbul','SBdisk','SBbul'], 25)

    df = mm.merge(meta[['Galaxy', 'L36', 'MHI', 'Reff']], on="Galaxy", how="inner")
    df = df[df["Galaxy"].isin(OUTLIERS)].dropna().copy()

    print(f"\n🚀 EXECUTING HIGH-MASS OUTLIER AUDIT")
    print("="*80)
    print(f"{'Galaxy':<12} | {'Ups':<6} | {'RMSE':<10} | {'Outer Res':<12} | {'Status'}")
    print("-" * 80)

    for gal in OUTLIERS:
        group = df[df["Galaxy"] == gal].copy()
        if group.empty: continue
        
        R, Vobs, errV = group["R"].values, group["Vobs"].values, group["errV"].values
        Vgas, Vdisk, Vbul = group["Vgas"].values, group["Vdisk"].values, group["Vbul"].values
        Reff_val, L36_val, MHI_val = group["Reff"].iloc[0], group["L36"].iloc[0], group["MHI"].iloc[0]
        
        weights = 1.0 / (errV**2 + SIGMA_INT**2 + SIGMA_STRUCT**2)
        best_chi2, best_pred, best_ups = np.inf, None, 0.5
        
        for ups in UPS_GRID:
            v_bar = np.sqrt(np.maximum(Vgas**2 + ups*Vdisk**2 + ups*Vbul**2, 1e-12))
            reff_m = max(float(Reff_val), 0.1) * C.KPC_TO_M
            sigma_si = (float(L36_val) * 1e9 * 1.988e30) / (np.pi * reff_m**2)
            
            vp = np.array([engine.predict_v_sparc(v_bar[k], sigma_si, Reff_val, L36_val, MHI_val, R_kpc=R[k]) for k in range(len(R))])
            
            chi = np.sum(weights * (Vobs - vp)**2)
            if chi < best_chi2:
                best_chi2, best_pred, best_ups = chi, vp, ups

        rmse = np.sqrt(np.sum(weights * (Vobs - best_pred)**2) / np.sum(weights))
        outer_res = Vobs[-1] - best_pred[-1]
        status = "✅ PASS" if abs(outer_res) < 60 else "⚠️ DEFICIT"
        
        print(f"{gal:<12} | {best_ups:<6.2f} | {rmse:<10.2f} | {outer_res:<+12.2f} | {status}")

    print("="*80)