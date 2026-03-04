#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - COMA CLUSTER INVARIANT INJECTION V.6.0
# Author: Miguel Navarro
# Logic: Viscous Transition in Isolated High-Mass Clusters.
# Protocol: Engine-Native Invariant Probe (Unified with EngineV6 Logic).
# =============================================================================

import numpy as np
import sys
import os

# Ensure the core is in the path
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, '../'))
sys.path.insert(0, PROJECT_ROOT)

try:
    from network_g_core import engineV6 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: Engine and Constants linked for Coma Probe.")
except ImportError:
    print("❌ ERROR: Could not find network_g_core. Ensure script is in project root.")
    sys.exit(1)

def run_isolated_coma_probe():
    print("══════════════════════════════════════════════════════════════════════════════════════════")
    print("🧪 ISOLATED COMA PROBE (ENGINE-NATIVE INJECTION)")
    print("══════════════════════════════════════════════════════════════════════════════════════════")

    # --- 1. Inputs (Locked Invariants for Coma Cluster) ---
    v_bar = 567.37
    logI_acc  = 0.038
    logI_comp = 3.92
    z_coma = 0.023
    v_target = 1441.12  # Observed velocity benchmark

    # --- 2. Call Engine Native Logic ---
    # Using predict_from_invariants ensures all Shatter Wall and 
    # Viscous Gain logic is handled automatically.
    v_pred, T = engine.predict_from_invariants(
        vbar_kms=v_bar, 
        logI_acc=logI_acc, 
        logI_comp=logI_comp, 
        redshift=z_coma, 
        return_trace=True
    )

    # --- 3. Statistical Analysis ---
    err = abs(v_pred - v_target) / v_target * 100.0
    boost = v_pred / v_bar

    # --- 4. Detailed Trace Reporting ---
    print(f"[INPUT]  v_bar:        {v_bar:.2f} km/s")
    print(f"[INPUT]  z_coma:       {z_coma:.3f}")
    print(f"[TRACE]  Visc Gain:    {T.get('visc_gain', 0):.6f}")
    print(f"[TRACE]  Eta Eff:      {T.get('eta_eff', 0):.6f}")
    print(f"[TRACE]  Shatter Supp: {T.get('shatter_suppress', 0):.6f}")
    print(f"[RESULT] Predicted V:  {v_pred:.2f} km/s")
    print(f"[RESULT] Target V:     {v_target:.2f} km/s")
    print(f"[RESULT] Engine Boost: {boost:.3f}")
    
    print("-" * 90)
    if err < 5.0:
        print(f"✅ PASS: Coma alignment within tolerance ({err:.2f}% deviation)")
    else:
        print(f"❌ FAIL: Coma deviation too high ({err:.2f}% deviation)")
    print("══════════════════════════════════════════════════════════════════════════════════════════")

if __name__ == '__main__':
    run_isolated_coma_probe()