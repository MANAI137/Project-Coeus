#!/usr/bin/env python3
# # =============================================================================
# PROJECT COEUS - CORNERSTONE SPARC AUDIT V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Vacuum Stiffening reclaiming R² = 0.956 / RMSE = 18.25 km/s.
# Protocol: 175 Galaxy Full-Sample Systematic Sweep.
# ================================
import numpy as np
import pandas as pd
import io, os, sys, time

# 1. SYSTEM LINKING
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
try:
    from network_g_core import engineV5 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: EngineV5 linked with V15.8 Statistics Bridge.")
except ImportError:
    print("❌ ERROR: Physics engine missing.")
    sys.exit(1)

DOI = "10.5281/zenodo.18641401"
ORCID = "0009-0009-5600-7985"

# CONFIG: THE HIGH-GROUND GATES
SIGMA_STRUCT, SIGMA_INT = 11.5, 3.0
INC_MIN_DEG, Q_MAX = 30.0, 2
UPS_GRID = np.linspace(0.4, 1.0, 25) 

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_CSV = os.path.join(DATA_DIR, "sparc_rotation_analyzed.csv")

def load_sparc_mrt(path, columns, skip_lines):
    if not os.path.exists(path):
        print(f"❌ ERROR: File not found at {path}")
        return pd.DataFrame(columns=columns)
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln for ln in f.readlines()[skip_lines:] if ln.strip() and not ln.strip().startswith(('#', '|'))]
    return pd.read_csv(io.StringIO("".join(lines)), sep=r"\s+", names=columns, engine="python", na_values=["----", "...", "NaN"])

if __name__ == "__main__":
    start_time = time.time()
    
    # 1. LOAD DATA
    meta_path = os.path.join(DATA_DIR, "sparc_metadata.csv")
    mm_path = os.path.join(DATA_DIR, "sparc_mass_models.csv")
    
    meta = load_sparc_mrt(meta_path, ['Galaxy','T','D','e_D','f_D','Inc','e_Inc','L36','e_L36','Reff','SBeff','Rdisk','SBdisk','MHI','RHI','Vflat','e_Vflat','Q','Ref'], 34)
    mm = load_sparc_mrt(mm_path, ['Galaxy','D_mm','R','Vobs','errV','Vgas','Vdisk','Vbul','SBdisk','SBbul'], 25)

    if meta.empty or mm.empty:
        sys.exit(1)

    # 2. SANITIZE AND FILTER (STRICT 149 GALAXY GATE)
    for col in ["Inc", "Q", "L36", "MHI", "Reff", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "R"]:
        if col in meta.columns: meta[col] = pd.to_numeric(meta[col], errors='coerce')
        if col in mm.columns: mm[col] = pd.to_numeric(mm[col], errors='coerce')

    df = mm.merge(meta[['Galaxy', 'Inc', 'Q', 'L36', 'MHI', 'Reff']], on="Galaxy", how="inner")
    
    # Apply filters: Inclination > 30, Quality <= 2, and non-null critical values
    df = df[(df["Q"] <= Q_MAX) & (df["Inc"] > INC_MIN_DEG)].dropna(subset=["Vobs", "Vgas", "R", "L36", "MHI", "Reff"]).copy()
    df["Vbul"] = df["Vbul"].fillna(0.0)

    unique_galaxies = df["Galaxy"].unique()
    print(f"🚀 EXECUTING IPP SPARC TRACE [Found {len(unique_galaxies)} Filtered Galaxies]")
    print("-" * 72)
    
    all_p, all_t, all_w, all_vbar = [], [], [], []
    analyzed_data_collection = []
    
    for gal, group in df.groupby("Galaxy"):
        group = group.copy()
        R, Vobs, errV = group["R"].values, group["Vobs"].values, group["errV"].values
        Vgas, Vdisk, Vbul = group["Vgas"].values, group["Vdisk"].values, group["Vbul"].values
        
        # Structural Anchors
        Reff_val = group["Reff"].iloc[0]
        L36_val = group["L36"].iloc[0]
        MHI_val = group["MHI"].iloc[0]
        
        # Statistical weights matching weighted R2 calculation
        weights = 1.0 / (errV**2 + SIGMA_INT**2 + SIGMA_STRUCT**2)
        best_chi2, best_pred, best_ups, best_vbar = np.inf, None, 0.52, None
        
        # --- Inside the Galaxy Loop ---
        for ups in UPS_GRID:
            # Calculate combined Baryonic velocity for current Ups
            v_bar = np.sqrt(np.maximum(Vgas**2 + ups*Vdisk**2 + ups*Vbul**2, 1e-12))
            reff_m = max(float(Reff_val), 0.1) * C.KPC_TO_M
            sigma_si = (float(L36_val) * 1e9 * 1.988e30) / (np.pi * reff_m**2)
            
            vp = np.array([
                engine.predict_v_sparc(
                    v_bar[k], 
                    sigma_si, 
                    Reff_val, 
                    L36_val, 
                    MHI_val,
                    R_kpc=R[k] 
                ) for k in range(len(R))
            ])
            
            chi = np.sum(weights * (Vobs - vp)**2)
            if chi < best_chi2:
                best_chi2, best_pred, best_ups, best_vbar = chi, vp, ups, v_bar

        if best_pred is not None:
            all_p.extend(best_pred); all_t.extend(Vobs); all_w.extend(weights); all_vbar.extend(best_vbar)
            group["Vpred"] = best_pred
            group["Vbar"] = best_vbar
            group["Residual"] = Vobs - best_pred
            group["Weight"] = weights 
            group["BestUps"] = best_ups
            analyzed_data_collection.append(group)
            
        print(f"✔️  Processed: {gal:<20}", end="\r")

    # 3. STATS AND OUTPUT
    y_p, y_t, w = np.array(all_p), np.array(all_t), np.array(all_w)
    
    # R-Squared (Weighted)
    ss_res = np.sum(w * (y_t - y_p)**2)
    y_mean = np.average(y_t, weights=w)
    ss_tot = np.sum(w * (y_t - y_mean)**2)
    r2 = 1 - (ss_res / ss_tot)
    
    # RMSE (km/s) - The target check for 17.83 km/s
    rmse_unweighted = np.sqrt(np.mean((y_t - y_p)**2))
    rmse_weighted = np.sqrt(np.sum(w * (y_t - y_p)**2) / np.sum(w))
    
    # 4. CSV EXPORT (STRICTLY FILTERED OUTPUT)
    if analyzed_data_collection:
        analyzed_df = pd.concat(analyzed_data_collection)
        analyzed_df = analyzed_df.rename(columns={'errV': 'eVobs'})
        
        if not os.path.exists(DATA_DIR):
            os.makedirs(DATA_DIR)
            
        analyzed_df.to_csv(OUTPUT_CSV, index=False)
        print(f"\n✅ DATA PERSISTED ({len(analyzed_data_collection)} Galaxies): {OUTPUT_CSV}")
    
    print("\n" + "="*72)
    print(f"FINAL GLOBAL WEIGHTED R² SCORE: {r2:.4f}")
    print(f"GLOBAL RMSE (Unweighted):      {rmse_unweighted:.2f} km/s")
    print(f"GLOBAL RMSE (Weighted):        {rmse_weighted:.2f} km/s")
    print(f"TOTAL POINTS PROCESSED:        {len(y_p)}")
    print(f"COMPLETION TIME:               {time.time()-start_time:.1f}s")
    print("="*72)