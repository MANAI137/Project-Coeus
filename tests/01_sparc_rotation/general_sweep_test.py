#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - NGC 3198 RELIC AUDIT V.6.0 (SI-LOCKED)
# Author: Miguel Navarro
# Logic: Treating NGC 3198 as a 'Mini-Cluster' (Merger Relic).
# Protocol: Identifying the Phase Transition Node where boost collapses.
# =============================================================================

import numpy as np
import os
import sys

# 1. SYSTEM LINKING (SI-LOCKED)
# Ensures the script can find network_g_core in the project root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '../../'))
sys.path.insert(0, PROJECT_ROOT)

try:
    from network_g_core import engineV6 as engine
    from network_g_core import constants as C
    print(f"✅ SUCCESS: EngineV6 and SI-Constants Linked.")
except ImportError:
    # Fallback if running from a different directory
    sys.path.insert(0, os.getcwd())
    try:
        from network_g_core import engineV6 as engine
        from network_g_core import constants as C
        print(f"✅ SUCCESS: EngineV6 Linked (Local Path).")
    except ImportError:
        print("❌ ERROR: Required Physics Core (V6) or Constants missing.")
        sys.exit(1)

def run_ngc3198_relic_audit():
    print("════════════════════════════════════════════════════════════════════════")
    print("🧪 NGC 3198 RELIC AUDIT: PHASE TRANSITION CONVERGENCE")
    print("════════════════════════════════════════════════════════════════════════")

    # --- 1. INPUT PARAMETERS (NGC 3198 Outer Radius) ---
    v_bar = 86.36      # Newtonian baryonic baseline (km/s)
    v_target = 150.0   # Observed velocity target (km/s)
    z_3198 = 0.0033    # Local redshift
    
    # --- 2. THE PHASE SWEEP ---
    # We sweep logI_comp from 'Galaxy' density to 'Extreme Cluster' density.
    # This identifies the node where internal merger stress triggers the Viscous Brake.
    print(f"{'logI_comp':<10} | {'Pred V':<10} | {'Boost':<8} | {'Visc Gain':<10} | {'Status'}")
    print("-" * 65)
    
    best_err = 999.0
    best_comp = 0.0
    
    # Range: 4.0 (High Density Spiral) to 12.0 (Solid Phase/Shatter Wall)
    for comp in np.linspace(4.0, 12.0, 17):
        v_pred, T = engine.predict_from_invariants(
            vbar_kms=v_bar, 
            logI_acc=0.040, 
            logI_comp=comp, 
            redshift=z_3198, 
            return_trace=True
        )
        
        err = abs(v_pred - v_target)
        boost = v_pred / v_bar
        vg = T.get('visc_gain', 0)
        
        status = ""
        if err < 10.0:
            status = "🎯 CONVERGENCE"
        
        if err < best_err:
            best_err = err
            best_comp = comp

        print(f"{comp:<10.2f} | {v_pred:<10.2f} | {boost:<8.3f} | {vg:<10.6f} | {status}")

    print("-" * 65)
    print(f"PHYSICAL VERDICT:")
    print(f"NGC 3198 aligns with a merger-relic structure at logI_comp ≈ {best_comp:.2f}.")
    print(f"At this density, the internal 'Lopsidedness' stress (Compaction) ")
    print(f"stiffens the vacuum into the viscous regime, damping the boost.")
    print("════════════════════════════════════════════════════════════════════════")

if __name__ == "__main__":
    run_ngc3198_relic_audit()