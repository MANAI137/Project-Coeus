import numpy as np
import os
from network_g_core import engineV6 as engine

def ngc5055_phase_split_resolution():
    # --- PHYSICAL DATA FOR NGC 5055 ---
    # R_25 (Optical Limit) is approx 11.6 kpc
    # Warp starts at 11.6 kpc and extends to 60+ kpc
    radii = [5.0, 10.0, 15.0, 25.0, 40.0, 55.0]
    v_obs = [190.0, 206.0, 185.0, 195.0, 205.0, 210.0] # Approx Sunflower curve
    v_bar = [150.0, 155.0, 110.0, 95.0, 85.0, 75.0]   # Newtonian Baseline

    print(f"--- NGC 5055 PHASE-SPLIT RESOLUTION ---")
    print(f"{'Radius':<8} | {'Phase':<12} | {'Node':<6} | {'Pred V':<8} | {'Error'}")
    print("-" * 60)

    for i in range(len(radii)):
        # 1. ZONE IDENTIFICATION
        if radii[i] < 11.6:
            # INNER CORE: "Gas Planet" Phase (Locked/Newtonian)
            # We use a very low logI_comp to represent the 'Locked' symmetry
            phase = "GAS PLANET"
            node = 1.0  
        else:
            # OUTER WARP: "Relic" Phase (Disturbed/Viscous)
            # We apply the 5.0 node we found in NGC 3198
            phase = "RELIC/WARP"
            node = 5.0

        v_pred = engine.predict_from_invariants(v_bar[i], 0.040, node, 0.0016)
        err = v_pred - v_obs[i]
        
        print(f"{radii[i]:<8.1f} | {phase:<12} | {node:<6.1f} | {v_pred:<8.2f} | {err:<+7.2f}")

if __name__ == "__main__":
    ngc5055_phase_split_resolution()