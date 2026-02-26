#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - PILLAR 3 BULLET SHATTER WALL V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Solid Phase Vacuum Saturation (1000 km/s Bump-Stop).
# Protocol: Engine Guardrail Phase-State Resolution (test_bullet.py).
# =============================================================================
import sys
import os

# 1. ENGINE LINKING
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
try:
    from network_g_core import engineV5 as engine
    print(f"✅ SUCCESS: Network G Engine (V5.9.0) linked.")
except ImportError:
    print("❌ ERROR: Physics engine missing.")
    sys.exit(1)

def run_shatter_test():
    print(f"\n🛡️  EXECUTING TEST 05: SHATTER WALL (V5.9.0)")
    print("-" * 65)

    # 3. DEFINE "SOLID PHASE" INVARIANTS (From Unified Harness)
    # High kinetic baseline (3000 km/s) to trip the wall
    vbar_in = 3000.0
    logI_acc = 0.0250
    logI_comp = 0.3010

    # 4. AUTHORITATIVE INJECTION
    v_pred, trace = engine.predict_from_invariants(
        vbar_kms=vbar_in,
        logI_acc=logI_acc,
        logI_comp=logI_comp,
        redshift=0.0,
        return_trace=True
    )

    # 5. VERIFICATION
    boost = v_pred / vbar_in
    target_boost = 1.00
    
    print(f"Baryonic Velocity: {vbar_in:>8.2f} km/s")
    print(f"Predicted Velocity: {v_pred:>8.2f} km/s")
    print(f"Engine Boost:       {boost:>8.3f}")
    print(f"Shatter Brake Vel:  {trace['shatter_brake_vel']:>8.6f}")
    print(f"Shatter Suppress:   {trace['shatter_suppress']:>8.6f}")
    print("-" * 65)

    # SUCCESS CRITERIA
    if abs(boost - target_boost) < 1e-4 and trace['shatter_brake_vel'] > 0.99:
        print("✅ PASS: Shatter Wall Guardrail Verified.")
    else:
        print("❌ FAIL: Boost not suppressed (Check V_SHATTER threshold).")
    print("=" * 65 + "\n")

if __name__ == "__main__":
    run_shatter_test()