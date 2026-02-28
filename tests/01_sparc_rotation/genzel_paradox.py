#!/usr/bin/env python3
# # =============================================================================
# PROJECT COEUS - PILLAR 2 HIGH-REDSHIFT DYNAMICS V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Expansion Damping & Baryon Dominance in the Early Universe (z~2).
# Protocol: Genzel Paradox Curve Simulation & Graphical Analysis.
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

# IPP V5 CONSTANTS
MODULUS_PERSISTENCE = 0.06113
I_THETA_LOCAL = 4.8e-10  # z = 0
I_THETA_Z2 = 1.2e-10     # z = 2 (Simulated Expansion Damping)
REFF = 5.0               # kpc

def get_ipp_curve(r, i_theta):
    # Newtonian Baseline (Simplified Exponential Disk)
    v_newton = 200 * (1 - np.exp(-r/2.0))
    
    # IPP Boost Logic
    # Saturated Phase (r <= Reff)
    boost_inner = MODULUS_PERSISTENCE * i_theta * (r * 3.086e19)
    # Persistent Phase (r > Reff) - Simplified for conceptual plot
    boost_outer = MODULUS_PERSISTENCE * i_theta * (REFF * 3.086e19) * (r/REFF)**-0.0605
    
    boost = np.where(r <= REFF, boost_inner, boost_outer)
    return np.sqrt(v_newton**2 + (boost / 1e6)) # Convert boost to km/s squared

r = np.linspace(0.1, 25, 100)
v_local = [get_ipp_curve(ri, I_THETA_LOCAL) for ri in r]
v_z2 = [get_ipp_curve(ri, I_THETA_Z2) for ri in r]
v_bar = 200 * (1 - np.exp(-r/2.0))

plt.figure(figsize=(10, 6))
plt.plot(r, v_local, 'b-', linewidth=2.5, label='Local Universe (z=0): Persistent Phase')
plt.plot(r, v_z2, 'm-', linewidth=2.5, label='Early Universe (z=2): Expansion Damped')
plt.plot(r, v_bar, 'k--', alpha=0.5, label='Baryonic Newtonian Limit')

plt.axvline(REFF, color='gray', linestyle=':', label='Effective Radius (Screening Gate)')
plt.title("The Genzel Paradox Resolution: Expansion-Driven Damping", fontsize=14, fontweight='bold')
plt.xlabel("Radius (kpc)")
plt.ylabel("Circular Velocity (km/s)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("IPP_Genzel_Paradox_V5.png", dpi=300)
plt.show()