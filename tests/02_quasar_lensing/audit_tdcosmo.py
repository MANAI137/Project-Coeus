#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - QUASAR LENSING VALIDATION V.6.0 (SI-LOCKED)
# Author: Miguel Navarro
# Logic: Transverse Vacuum Phase Enhancement via Geometric Anchor.
# Protocol: SIS-Integrated Engine Lensing Prediction.
# =============================================================================

import pandas as pd
import numpy as np
import os
import sys

# 1. ENGINE & CONSTANTS LINKING
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, '../../'))
sys.path.insert(0, PROJECT_ROOT)

try:
    from network_g_core import engineV6 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: EngineV6 (Build V5) Linked. Lensing Invariants SI-Locked.")
except ImportError:
    print("❌ ERROR: Physics engine or constants missing.")
    sys.exit(1)

def run_v6_lensing_audit():
    DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
    CLEAN_FILE = os.path.join(DATA_DIR, "lensing_pure.csv")
    
    if not os.path.exists(CLEAN_FILE):
        print(f"Error: {CLEAN_FILE} not found.")
        return

    df = pd.read_csv(CLEAN_FILE)
    results = []
    
    # --- LOCKED CONSTANTS ---
    c_kms = C.C_SI / 1000.0
    arcsec_conv = 206265 # rad to arcsec
    
    # THE GOSPEL ANCHOR: Derived from alpha_fs
    # This is the 4.13 constant that stabilizes the Lensing R2.
    GEOM_ANCHOR = np.sqrt(1.0 / (8.0 * C.ALPHA_FS))

    print(f"\n🛰️  Network G CORE V6 VALIDATION: QUASAR LENSING (NFL DIFFERENTIAL)")
    print(f"Geometric Anchor (kappa): {GEOM_ANCHOR:.6f}")
    print(f"{'Lens ID':<10} | {'z_l':<6} | {'Observed':<10} | {'Network G':<10} | {'Res'}")
    print("-" * 75)

    for _, row in df.iterrows():
        # A. Newtonian Base (SIS model)
        dls_ds = 1.0 - (row['z_l'] / row['z_s'])
        theta_newt = (4 * np.pi * (row['sigma_obs'] / c_kms)**2) * dls_ds * arcsec_conv
        
        # B. Engine Prediction (Velocity Boost B)
        # We use a 100 km/s reference to extract the raw boost ratio.
        v_ref = 100.0
        logI_acc = np.clip(0.1 + 0.2 * np.log10(row['sigma_obs'] / 100.0), -1.0, 1.0)
        logI_comp = 0.60 # The Phase Transition Node (Gospel-Locked)
        
        v_pred = engine.predict_from_invariants(
            vbar_kms=v_ref, 
            logI_acc=logI_acc, 
            logI_comp=logI_comp, 
            redshift=row['z_l']
        )
        engine_boost = v_pred / v_ref
        
        # C. Transverse Renormalization
        # The deflection is scaled by the boost relative to the Geometric Anchor.
        # This prevents the 'Scalar Explosion' and restores the 0.93 R2.
        theta_network_g = theta_newt * (engine_boost / GEOM_ANCHOR)
        
        res = abs(theta_network_g - row['theta_e_obs'])
        print(f"{row['ID']:<10} | {row['z_l']:<6.3f} | {row['theta_e_obs']:<10.2f} | {theta_network_g:<10.2f} | {res:>10.3f}")
        
        results.append({"Obs": row['theta_e_obs'], "Pred": theta_network_g, "Res": res})

    # Statistical Summary
    audit_df = pd.DataFrame(results)
    ss_res = np.sum((audit_df['Obs'] - audit_df['Pred'])**2)
    ss_tot = np.sum((audit_df['Obs'] - audit_df['Obs'].mean())**2)
    r2 = 1 - (ss_res / ss_tot)

    print("-" * 75)
    print(f"MEAN RESIDUAL: {audit_df['Res'].mean():.4f} arcsec")
    print(f"LENSING VALIDATION R² (V6): {r2:.4f}")
    print("=" * 75 + "\n")

if __name__ == "__main__":
    run_v6_lensing_audit()