#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - PILLAR 3 GHOST DISCRIMINATION V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Invariant A/B Structural Divergence (Gas vs. Galaxies).
# Protocol: 2.33x Boost Ratio Verification (test_bullet_ghost.py).
# =============================================================================

import sys
import os

# 1. ENGINE LINKING
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
try:
    from network_g_core import engineV5 as engine
    print(f"✅ SUCCESS: Network G Engine (V5.9.0) linked.")
except ImportError:
    print("❌ ERROR: Physics engine (network_g_core/engineV5.py) missing.")
    sys.exit(1)

def run_ghost_audit():
    print(f"\n🛰️  EXECUTING TEST 06: BULLET CLUSTER GHOST (V5.9.0)")
    print("-" * 65)

    # 3. DEFINE INVARIANTS (From Unified Harness V11)
    # Point A: Baryonic Gas (High density, low network connectivity)
    # Point B: Galaxy Center (Lower density, high structural network stress)
    test_points = {
        "Point A (Gas)":      {"logI_acc": 0.0200, "logI_comp": 3.5000, "vbar": 567.37},
        "Point B (Galaxies)": {"logI_acc": 0.0500, "logI_comp": 0.5000, "vbar": 567.37}
    }

    results = {}

    for name, data in test_points.items():
        # Using authoritative invariant injection path
        v_pred, trace = engine.predict_from_invariants(
            vbar_kms=data["vbar"],
            logI_acc=data["logI_acc"],
            logI_comp=data["logI_comp"],
            redshift=0.0,
            return_trace=True
        )
        
        boost = v_pred / data["vbar"]
        results[name] = boost
        
        print(f"{name:<20} | Boost: {boost:.3f} | acc: {trace['logI_acc']:.4f} | comp: {trace['logI_comp']:.4f}")

    # 4. RATIO VERIFICATION
    ratio = results["Point B (Galaxies)"] / results["Point A (Gas)"]
    target_ratio = 2.334
    fidelity = (1.0 - abs(ratio - target_ratio) / target_ratio) * 100

    print("-" * 65)
    print(f"RESULTING BOOST RATIO (B/A): {ratio:.3f}")
    print(f"TARGET RATIO (Unified V11): {target_ratio:.3f}")
    print(f"STRUCTURAL FIDELITY:        {fidelity:.2f}%")
    
    if fidelity > 99.0:
        print("\n✅ PASS: Network G Discrimination Verified.")
    else:
        print("\n❌ FAIL: Disconnect in Structural Invariants.")
    print("=" * 65 + "\n")

if __name__ == "__main__":
    run_ghost_audit()