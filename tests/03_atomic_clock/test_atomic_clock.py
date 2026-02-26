#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - ATOMIC CLOCK NULL TEST V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Engine-Integrated Shielding (rho_u -> S_shield).
# Protocol: Galileo E05 SP3 Eccentricity validation against GREAT (2018).
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd

# 1. SYSTEM LINKING
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
try:
    from network_g_core import engineV5 as engine
    from network_g_core import constants
    print(f"✅ SUCCESS: Network G Engine (V5.9.0) linked.")
except ImportError:
    print("❌ ERROR: Physics engine missing.")
    sys.exit(1)


# 3. SOLAR SYSTEM PARAMETERS
GM_EARTH = 3.986004418e14  # m^3/s^2
C_LIGHT  = 299792458.0     # m/s

def parse_sp3_fixed(filename, sat_id="E05"):
    """Fixed-width parser for SP3 GPS position records."""
    data = []
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Missing orbital data: {filename}")
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("P") and sat_id in line[:8]:
                try:
                    x, y, z = float(line[4:18])*1000, float(line[18:32])*1000, float(line[32:46])*1000
                    r = np.sqrt(x**2 + y**2 + z**2)
                    if np.isfinite(r) and r > 0:
                        data.append((r, np.sqrt(GM_EARTH / r)))
                except (ValueError, IndexError): continue
    return np.asarray(data)

def run_clock_audit(xi_unscreened=1.09, rho_crit=1e-18, p=8.0, great_frac=2.49e-5):
    print(f"\n🛰️  EXECUTING TEST 04: ATOMIC CLOCK NULL-TEST (V5.9.0)")
    print("-" * 88)

    sp3_path = os.path.join(os.path.dirname(__file__), "data", "esoc20446.sp3")
    orbit_data = parse_sp3_fixed(sp3_path)
    r, v = orbit_data[:, 0], orbit_data[:, 1]
    r_mean = np.mean(r)

    # 4. GR BASELINE & KINETIC TERMS
    phi_local = GM_EARTH / (r * C_LIGHT**2)
    phi_ref   = GM_EARTH / (r_mean * C_LIGHT**2)
    kinetic   = (v**2) / (2.0 * C_LIGHT**2)
    y_gr = -(phi_local - phi_ref) - kinetic

    # 5. ENGINE SHIELDING CALCULATION
    # Local field energy density: rho_u = g^2 / (8*pi*G*c^2)
    g_local = GM_EARTH / (r**2)
    rho_u = (g_local**2) / (8.0 * np.pi * constants.G_SI * (C_LIGHT**2))
    
    # Authoritative Shielding Law
    S = 1.0 / (1.0 + (rho_u / rho_crit) ** p)

    # 6. NUMERICAL STABILITY FIX: Compute shift directly
    # This prevents float64 rounding from zeroing out the tiny Network G signal
    xi_shift = (xi_unscreened - 1.0) * S
    xi_eff   = 1.0 + xi_shift 
    y_network_g = -(phi_local - phi_ref) * xi_eff - kinetic

    # 7. VERDICT VERIFICATION
    max_abs_err = np.max(np.abs(y_network_g - y_gr))
    threshold = (np.max(y_gr) - np.min(y_gr)) * great_frac

    print(f"{'Max |Network G - GR| deviation':<60} | {max_abs_err:<20.12e}")
    print(f"{'GREAT (2018) detectability limit':<60} | {threshold:<20.6e}")
    print("-" * 88)

    if max_abs_err < threshold:
        print("VERDICT: PASS ✅ (Network G is hidden by Shielding Law)")
    else:
        print("VERDICT: FAIL ❌ (Experimental violation detected)")
    print("=" * 88 + "\n")

if __name__ == "__main__":
    run_clock_audit()