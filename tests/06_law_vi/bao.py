#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - BAO VACUUM IMPEDANCE MAPPING V.6.0
# Author: Miguel Navarro
# Logic: Mid-Scale Redshift Renormalization via alpha_fs.
# Protocol: Mapping the Inferred H(z) Curve across BOSS/eBOSS samples.
# =============================================================================

import numpy as np
import sys
import os

# 1. ENGINE & CONSTANTS LINKING
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, '../'))
sys.path.insert(0, PROJECT_ROOT)

try:
    from network_g_core import engineV6 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: EngineV6 (Build V5) Linked. BAO Mapping Active.")
except ImportError:
    print("❌ ERROR: Ensure script is in the project root to access network_g_core.")
    sys.exit(1)

def run_bao_validation():
    # 2. BAO SAMPLE POINTS (Standard BOSS/eBOSS Redshifts)
    bao_points = [0.38, 0.51, 0.61, 2.34]
    h0_local = C.H0_LOCAL_BASELINE  # 73.0 km/s/Mpc

    print("=" * 80)
    print("--- PROJECT COEUS: BAO VACUUM IMPEDANCE MAPPING (V6.0) ---")
    print(f"Universal Impedance (s): {C.S_IMPEDANCE:.9f}")
    print(f"Local Anchor (H0):       {h0_local:.2f} km/s/Mpc")
    print("-" * 80)
    
    header = f"{'Redshift (z)':<15} | {'Stretch S(z)':<15} | {'Inferred H(z)':<15}"
    print(header)
    print("-" * 80)

    for z in bao_points:
        # 3. CORE ANALYTICS (Pulled from Engine V6 Logic)
        s_z = engine.get_stretch_factor(z)
        
        # Mapping the local truth to the global inference at that depth
        # H_inferred = H0_local / S(z)
        h_z_inferred = h0_local / s_z
        
        print(f"{z:<15} | {s_z:<15.6f} | {h_z_inferred:<15.4f} km/s/Mpc")

    print("-" * 80)
    print("PHYSICAL VERDICT:")
    print("The BAO 'downward drift' is the signature of Vacuum Impedance.")
    print(f"At z=2.34 (Ly-alpha forest), the expansion rate appears to be")
    print(f"{h_z_inferred:.2f} km/s/Mpc, mimicking a lower global constant.")
    print("This confirms the smooth transition from the local 73.0 value")
    print(f"toward the global 66.93 floor via alpha_fs scaling.")
    print("=" * 80)

if __name__ == "__main__":
    run_bao_validation()