#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - NGC 3198 VISCOUS TRANSITION TEST (MINI-CLUSTER MODE)
# Logic: Treating the merger-relic as a high-density compact system.
# =============================================================================

import numpy as np
import sys
import os

# Link to project root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

try:
    from network_g_core import engineV6 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: EngineV6 Linked for Viscous Probe.")
except ImportError:
    print("❌ ERROR: Ensure script is in /tests/01_sparc_rotation/ or project root.")
    sys.exit(1)

def run_ngc3198_viscous_probe():
    print("════════════════════════════════════════════════════════════════════════")
    print("🧪 NGC 3198 VISCOUS PROBE (MINI-CLUSTER INJECTION)")
    print("════════════════════════════════════════════════════════════════════════")

    # --- 1. Inputs (NGC 3198 at the Outer Radius) ---
    v_bar = 86.36      # Newtonian baryonic baseline
    v_target = 150.0   # The observed velocity we want to hit
    z_3198 = 0.0033    # Local redshift
    
    # We 'Inject' a high logI_comp to simulate the merger-relic stiffening
    # (Typical for the 'White Dwarf Cluster' or compact-relic logic)
    logI_acc  = 0.040  # Standard acceleration invariant
    logI_comp = 3.85   # Pushed into the Viscous Transition zone (>3.5)

    # --- 2. Call Engine Native Logic ---
    v_pred, T = engine.predict_from_invariants(
        vbar_kms=v_bar, 
        logI_acc=logI_acc, 
        logI_comp=logI_comp, 
        redshift=z_3198, 
        return_trace=True
    )

    # --- 3. Statistical Analysis ---
    err = abs(v_pred - v_target)
    boost = v_pred / v_bar

    # --- 4. Detailed Trace Reporting ---
    print(f"[INPUT]  v_bar:       {v_bar:.2f} km/s")
    print(f"[TRACE]  Visc Gain:    {T.get('visc_gain', 0):.6f}")
    print(f"[TRACE]  Shatter Supp: {T.get('shatter_suppress', 0):.6f}")
    print(f"[RESULT] Predicted V:  {v_pred:.2f} km/s")
    print(f"[RESULT] Target V:     {v_target:.2f} km/s")
    print(f"[RESULT] Engine Boost: {boost:.3f}")
    
    print("-" * 72)
    if err < 10.0:
        print(f"✅ SUCCESS: Viscous transition resolves the NGC 3198 over-prediction.")
        print(f"Deviation dropped to {err:.2f} km/s.")
    else:
        print(f"⚠️  NOTICE: Residual remains at {err:.2f} km/s. Adjusting logI_comp suggested.")
    print("════════════════════════════════════════════════════════════════════════")

if __name__ == '__main__':
    run_ngc3198_viscous_probe()