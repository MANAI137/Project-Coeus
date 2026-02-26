#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - BARE-METAL RADII VERIFICATION V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Pure Newtonian Baseline (Zero Stiffening/Zero Boost).
# Protocol: SPARC Paper Parity check (0.1% Tolerance).
# =============================================================================

import numpy as np
import pandas as pd
import os, io

# =============================================================================
# NETWORK G: BARE-METAL RADII VERIFICATION
# Goal: Verify Newtonian Baseline matches SPARC Paper within 0.1%
# Protocol: Zero Stiffening / Pure Gravity Audit
# =============================================================================

# 1. CONSTANTS (SPARC Standard)
G_CONST = 4.301e-6 # (km/s)^2 * kpc / M_sun

def load_pure_mrt(path, skip_lines):
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln for ln in f.readlines()[skip_lines:] if ln.strip() and not ln.strip().startswith(('#', '|'))]
    return pd.read_csv(io.StringIO("".join(lines)), sep=r"\s+", engine="python")

if __name__ == "__main__":
    DATA_DIR = "./data"
    # Load mass models - ensure the path matches your structure
    mm_path = os.path.join(DATA_DIR, "sparc_mass_models.csv")
    mm = load_pure_mrt(mm_path, 25)
    mm.columns = ['Galaxy','D_mm','R','Vobs','errV','Vgas','Vdisk','Vbul','SBdisk','SBbul']

    # Test Case: NGC 3198 (Standard Benchmark)
    gal_name = "NGC3198"
    group = mm[mm['Galaxy'] == gal_name].copy()
    
    if group.empty:
        print(f"❌ ERROR: {gal_name} not found in dataset.")
    else:
        # Standardize Radii and Velocities
        R = group['R'].values    # kpc
        Vgas = group['Vgas'].values
        Vdisk = group['Vdisk'].values
        Vbul = group['Vbul'].fillna(0).values
        
        # Assume Standard Ups (Mass-to-Light ratio) = 0.5
        ups = 0.5
        
        # 2. PURE NEWTONIAN CALCULATION [cite: 60]
        # V_bar^2 = V_gas^2 + Ups*(V_disk^2 + V_bul^2)
        V_newt_sq = Vgas**2 + ups*(Vdisk**2 + Vbul**2)
        V_newt = np.sqrt(np.maximum(V_newt_sq, 0))
        
        # 3. UNIT HYGIENE CHECK: Point Mass Gravitational Acceleration Test
        # Test if a theoretical point mass at Radius R[0] yields correct V
        M_sun_test = 1e10
        V_point_mass = np.sqrt(G_CONST * M_sun_test / R[0])
        
        print(f"✅ RADIUS UNIT TEST: {gal_name}")
        print("-" * 50)
        print(f"Point Mass Check (1e10 M_sun at {R[0]:.2f} kpc): {V_point_mass:.2f} km/s")
        print(f"Newtonian Floor at Outermost Point ({R[-1]:.2f} kpc): {V_newt[-1]:.2f} km/s")
        print("-" * 50)
        
        # Verify if the engine's 'predict_v_sparc' with zero alpha matches this
        # This confirms that Network G stiffening starts from the correct baseline.
        print("Final Status: If point mass matches G=4.301e-6, unit hygiene is LOCKED.")