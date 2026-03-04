#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - NGC 3198 HIERARCHY AUDIT V.6.0 (SI-LOCKED)
# Protocol: 25% Anchor / 75% Blind Prediction (Relic Mode: logI_comp = 5.0)
# =============================================================================
import numpy as np
import pandas as pd
import io, os, sys

# 1. ENGINE LINKING
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '../../'))
sys.path.insert(0, PROJECT_ROOT)

try:
    from network_g_core import engineV6 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: EngineV6 Linked. Initiating Hierarchy Protocol.")
except ImportError:
    print("❌ ERROR: Physics engine or Constants missing.")
    sys.exit(1)

# CONFIG
SIGMA_STRUCT, SIGMA_INT = 11.5, 3.0
TARGET = "NGC3198"
RELIC_COMPACTION = 5.0 

def load_sparc_mrt(path, columns, skip_lines):
    if not os.path.exists(path): return pd.DataFrame()
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln for ln in f.readlines()[skip_lines:] if ln.strip() and not ln.strip().startswith(('#', '|'))]
    return pd.read_csv(io.StringIO("".join(lines)), sep=r"\s+", names=columns, engine="python", na_values=["----", "...", "NaN"])

def run_hierarchy_audit():
    DATA_DIR = os.path.join(SCRIPT_DIR, "data")
    META_PATH = os.path.join(DATA_DIR, "sparc_metadata.csv")
    MM_PATH = os.path.join(DATA_DIR, "sparc_mass_models.csv")

    meta = load_sparc_mrt(META_PATH, ['Galaxy','T','D','e_D','f_D','Inc','e_Inc','L36','e_L36','Reff','SBeff','Rdisk','SBdisk','MHI','RHI','Vflat','e_Vflat','Q','Ref'], 34)
    mm = load_sparc_mrt(MM_PATH, ['Galaxy','D_mm','R','Vobs','errV','Vgas','Vdisk','Vbul','SBdisk','SBbul'], 25)

    df = mm.merge(meta[['Galaxy', 'L36', 'MHI', 'Reff']], on="Galaxy", how="inner")
    group = df[df["Galaxy"] == TARGET].copy()

    if group.empty:
        print(f"❌ ERROR: Data for {TARGET} not found.")
        return

    # 2. HIERARCHY PROTOCOL: OUTSIDE-IN SORT (Descending Radius)
    group = group.sort_values("R", ascending=False)
    R, Vobs, errV = group["R"].values, group["Vobs"].values, group["errV"].values
    Vgas, Vdisk, Vbul = group["Vgas"].values, group["Vdisk"].values, group["Vbul"].values
    
    # Split: 25% Anchor Points (Outer)
    n_split = max(int(len(R) * 0.25), 2)
    w_pts = 1.0 / (errV**2 + SIGMA_INT**2 + SIGMA_STRUCT**2)

    print(f"\n🛰️  HIERARCHY AUDIT: {TARGET} (Relic Node: {RELIC_COMPACTION})")
    print(f"Anchor: {n_split} outer points | Blind: {len(R)-n_split} inner points")
    print("-" * 75)

    # 3. PREDICTION LOOP
    valid_p, valid_t, valid_w = [], [], []
    v_bar_all = np.sqrt(np.maximum(Vgas**2 + 1.0*Vdisk**2 + 1.0*Vbul**2, 1e-12))

    for j in range(len(R)):
        # Calculate logI_acc based on Newtonian baseline
        logI_acc = np.log10(np.maximum((v_bar_all[j]**2) / (R[j] * C.KPC_TO_M), 1e-15))
        
        # Calling engine without unpacking a second value to avoid TypeError
        vp = engine.predict_from_invariants(
            vbar_kms=v_bar_all[j],
            logI_acc=logI_acc,
            logI_comp=RELIC_COMPACTION, 
            redshift=0.0033
        )
        
        # Statistic collection for the 'Blind' zone (the inner 75%)
        if j >= n_split:
            valid_p.append(vp)
            valid_t.append(Vobs[j])
            valid_w.append(w_pts[j])

    # 4. RESULTS REPORTING
    y_p, y_t, w = np.array(valid_p), np.array(valid_t), np.array(valid_w)
    avg_res = np.mean(np.abs(y_t - y_p))
    rmse = np.sqrt(np.mean((y_t - y_p)**2))
    
    print(f"BLIND PREDICTION RESULTS (Inner 75%):")
    print(f"Mean Absolute Residual: {avg_res:.2f} km/s")
    print(f"RMSE (Blind Zone):      {rmse:.2f} km/s")
    print("-" * 75)
    
    if avg_res < 20.0:
        print("VERDICT: Hierarchy Logic VALIDATED.")
        print("The Relic Node predicts the blind zone with high fidelity.")
    else:
        print("VERDICT: Discrepancy detected. Reviewing Phase Transition bounds.")
    print("=" * 75 + "\n")

if __name__ == "__main__":
    run_hierarchy_audit()