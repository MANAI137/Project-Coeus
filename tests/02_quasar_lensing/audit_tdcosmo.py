#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - QUASAR LENSING AUDIT V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Vacuum Phase Enhancement in High-Redshift SIS Systems.
# Protocol: SIS-Integrated Engine Lensing Prediction (TDCOSMO Sync).
# =============================================================================

import pandas as pd
import numpy as np
import os
import sys
import io

# 1. SYSTEM LINKING
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
try:
    from network_g_core import engineV5 as engine
    from network_g_core import constants
    print(f"✅ SUCCESS: Network G Engine (V5.9.0) and Constants linked.")
except ImportError:
    print("❌ ERROR: Physics engine or constants missing.")
    sys.exit(1)

# 2. DEFINITIVE PROJECT IDENTIFIERS
DOI = "10.5281/zenodo.18641401"
ORCID = "0009-0009-5600-7985"

# 3. DIRECTORY CONFIG
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
CLEAN_FILE = os.path.join(DATA_DIR, "lensing_pure.csv")

def run_v5_audit():
    if not os.path.exists(CLEAN_FILE):
        print(f"Error: {CLEAN_FILE} not found.")
        return

    df = pd.read_csv(CLEAN_FILE)
    results = []
    
    # Constants from Gospel/SI
    c_kms = constants.C_SI / 1000.0
    arcsec_conv = 206265 # radians to arcseconds

    print(f"\n🛰️  Network G CORE V5 VALIDATION: QUASAR LENSING (NFL DIFFERENTIAL)")
    print(f"Ref DOI: {DOI} | Author: {ORCID}")
    print(f"{'Lens ID':<10} | {'z_l':<6} | {'Observed':<8} | {'Network G':<8} | {'Res'}")
    print("-" * 72)

    for _, row in df.iterrows():
        # A. Newtonian Base (SIS model)
        # Standard light bending: theta = 4 * pi * (sigma/c)^2 * (Dls/Ds)
        dls_ds = 1.0 - (row['z_l'] / row['z_s'])
        theta_newt = (4 * np.pi * (row['sigma_obs'] / c_kms)**2) * dls_ds * arcsec_conv
        
        # B. V5 Invariant Mapping
        # We extract the boost (B) from the engine using a reference baryonic velocity.
        v_ref = 100.0 
        # Quasar lenses are high-stress environments; logI_acc scales with dispersion.
        logI_acc = np.clip(0.1 + 0.2 * np.log10(row['sigma_obs'] / 100.0), -1.0, 1.0)
        # Use the Gospel Phase Transition Node (0.60 kg/m2)
        logI_comp = 0.60 
        
        v_pred = engine.predict_from_invariants(
            vbar_kms=v_ref, 
            logI_acc=logI_acc, 
            logI_comp=logI_comp, 
            redshift=row['z_l']
        )
        
        # C. Differential Enhancement
        # In V5, SPARC anchor is 4.13. We scale the SIS model by the engine's predicted boost.
        engine_boost = v_pred / v_ref
        
        # We normalize the boost to the expected gravitational enhancement.
        # This prevents the "Scalar Explosion" from the previous run.
        theta_network_g = theta_newt * (engine_boost / 4.13)
        
        res = abs(theta_network_g - row['theta_e_obs'])
        print(f"{row['ID']:<10} | {row['z_l']:<6.3f} | {row['theta_e_obs']:<8.2f} | {theta_network_g:<8.2f} | {res:>8.3f}")
        
        results.append({
            "ID": row['ID'],
            "Obs": row['theta_e_obs'], 
            "Pred": theta_network_g, 
            "Res": res
        })

    # Statistical Summary
    audit_df = pd.DataFrame(results)
    ss_res = np.sum((audit_df['Obs'] - audit_df['Pred'])**2)
    ss_tot = np.sum((audit_df['Obs'] - audit_df['Obs'].mean())**2)
    r2 = 1 - (ss_res / ss_tot)

    print("-" * 72)
    print(f"MEAN RESIDUAL: {audit_df['Res'].mean():.4f} arcsec")
    print(f"LENSING VALIDATION R² (V5): {r2:.4f}")
    print("=" * 72 + "\n")

if __name__ == "__main__":
    run_v5_audit()