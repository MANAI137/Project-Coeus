#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - SPARC GLOBAL DIAGNOSTICS V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Statistical Error Analysis & Master Correlation Mapping.
# Protocol: RMSE and Chi-Squared Validation across 4-Regime Invariants.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os
from scipy.optimize import minimize_scalar

# ---------------------------------------------------------
# IPP ENGINE V5 CONSTANTS (Gospel Parity)
# ---------------------------------------------------------
MODULUS_PERSISTENCE = 0.06113  # Alpha (Universal Modulus)
IMPEDANCE_SCALAR = 4.8e-10    # I_theta (Impedance Scalar)
GAMMA = -0.0605              # Impedance Exponent
G_SI = 6.67430e-11           # Gravitational Constant
MSUN_SI = 1.98847e30         # Solar Mass
KPC_TO_M = 3.085677581e19    # kpc to meters
SIGMA_REF = (100000.0**2) / (10.0 * KPC_TO_M) # Stress Reference (100km/s @ 10kpc)

# PATH LOGIC: Pin to script directory (/tests/01_sparc_rotation/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(SCRIPT_DIR, "data", "sparc_rotation_analyzed.csv")
OUT_DIR = os.path.join(SCRIPT_DIR, "visualizations")

# ---------------------------------------------------------
# CORE PHYSICS: INSIDE-OUT PREDICTION
# ---------------------------------------------------------
def compute_rscale_m(L36, MHI, ups_global):
    """Calculates the Metric Saturation Horizon (R_scale)."""
    Mbar_msun = 1.33 * MHI + ups_global * L36
    Mbar_kg = Mbar_msun * MSUN_SI
    return math.sqrt((G_SI * Mbar_kg) / (MODULUS_PERSISTENCE * IMPEDANCE_SCALAR))

def predict_v_ipp_kms(R_kpc, Reff_kpc, Vgas, Vdisk, Vbul, L36, MHI, ups):
    """IPP Engine V5: Inside-Out Velocity Prediction."""
    R_m = np.maximum(1e-6, R_kpc * KPC_TO_M)
    Reff_m = Reff_kpc * KPC_TO_M
    
    # Baryonic Potential (Newtonian)
    s = math.sqrt(max(0, ups))
    Vbar_ms = math.sqrt(max(0, Vgas**2 + (s*Vdisk)**2 + (s*Vbul)**2)) * 1e3
    
    # Saturation Horizon
    Rscale_m = compute_rscale_m(L36, MHI, ups)
    
    # Local Stress Profile
    Vdisk_scaled_ms = s * Vdisk * 1e3
    sigma_loc = (Vdisk_scaled_ms**2) / R_m
    sigma_ratio = max(1e-12, sigma_loc / SIGMA_REF)
    
    # Phase Transition Logic
    if R_m <= Reff_m:
        boost_v2 = MODULUS_PERSISTENCE * IMPEDANCE_SCALAR * R_m
    else:
        boost_v2 = (MODULUS_PERSISTENCE * IMPEDANCE_SCALAR * Rscale_m) * (sigma_ratio ** GAMMA)
        
    return math.sqrt(max(0, Vbar_ms**2 + boost_v2)) / 1e3

# ---------------------------------------------------------
# PLOT: MASTER CORRELATION (V_obs vs V_pred)
# ---------------------------------------------------------
def plot_master_correlation(df):
    """Generates the scattergram with Weighted R2 parity for Document B."""
    print("Generating Weighted Master Correlation Scattergram...")
    
    v_obs = df['Vobs'].values
    v_pred = df['Vpred'].values 
    weights = df['Weight'].values  # Imported from test_sparcV5 for mathematical parity
    
    n_galaxies = len(df['Galaxy'].unique())
    n_points = len(df)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Identity Line
    lims = [0, max(max(v_obs), max(v_pred)) + 20]
    ax.plot(lims, lims, color='gray', linestyle='--', alpha=0.8, label='Identity (1:1)', zorder=1)
    
    # Scatter points colored by Residual
    residuals = v_obs - v_pred
    sc = ax.scatter(v_obs, v_pred, c=residuals, cmap='RdBu_r', s=25, alpha=0.6, edgecolors='k', linewidth=0.5, zorder=2)
    
    # CALCULATE WEIGHTED R2 (Absolute Parity with test_sparcV5 logic)
    res_sq = weights * (v_obs - v_pred)**2
    tot_sq = weights * (v_obs - np.average(v_obs, weights=weights))**2
    r2_weighted = 1 - (np.sum(res_sq) / np.sum(tot_sq))
    rmse_weighted = np.sqrt(np.sum(res_sq) / np.sum(weights))
    
    # Stats box
    stats_text = f"$R^2_{{weighted}} = {r2_weighted:.4f}$\n" + rf"$RMSE_{{weighted}} = {rmse_weighted:.2f} \text{{ km/s}}$"
    
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, verticalalignment='top', 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel(r'$V_{obs}$ (km/s)', fontsize=12)
    ax.set_ylabel(r'$V_{pred}$ (km/s)', fontsize=12)
    ax.set_title(f'IPP Engine V5 Correlation ($N={n_galaxies}$ Galaxies, {n_points} Points)', fontsize=14, weight='bold')
    
    plt.colorbar(sc, label=r'Residual ($V_{obs} - V_{pred}$)')
    ax.grid(True, linestyle='--', alpha=0.3)
    
    out_path = os.path.join(OUT_DIR, 'IPP_Master_Correlation_N149.png')
    plt.savefig(out_path, dpi=300)
    print(f"Scattergram Saved -> {out_path}")

# ---------------------------------------------------------
# PLOT: ENERGY DENSITY AUDIT (Diagnostic)
# ---------------------------------------------------------
def plot_energy_density_audit(df, galaxy_name='F571-8'):
    """Diagnostic to audit Local vs Global Modulus alignment."""
    print(f"Auditing Energy Density Phase for {galaxy_name}...")
    g_data = df[df['Galaxy'] == galaxy_name].sort_values('R')
    if g_data.empty: return

    R_vals = g_data['R'].values
    Vobs_vals = g_data['Vobs'].values
    Vgas_vals = g_data['Vgas'].values
    Vdisk_vals = g_data['Vdisk'].values
    Vbul_vals = g_data['Vbul'].values
    Reff = g_data['Reff'].iloc[0]
    L36 = g_data['L36'].iloc[0]
    MHI = g_data['MHI'].iloc[0]
    
    implied_ups, implied_vpred = [], []
    for i in range(len(R_vals)):
        def objective(u):
            vp = predict_v_ipp_kms(R_vals[i], Reff, Vgas_vals[i], Vdisk_vals[i], Vbul_vals[i], L36, MHI, u)
            return (vp - Vobs_vals[i])**2
        res = minimize_scalar(objective, bounds=(0.01, 5.0), method='bounded')
        implied_ups.append(res.x)
        implied_vpred.append(predict_v_ipp_kms(R_vals[i], Reff, Vgas_vals[i], Vdisk_vals[i], Vbul_vals[i], L36, MHI, res.x))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    ax1.errorbar(R_vals, Vobs_vals, yerr=g_data['eVobs'], fmt='ko', label='Observed Data', zorder=5)
    ax1.plot(R_vals, implied_vpred, 'r-', linewidth=2.5, label=r'IPP Phase Trace', zorder=4)
    ax1.set_ylabel(r'Velocity (km/s)')
    ax1.set_title(f'{galaxy_name}: Phase Persistence Diagnostic', fontsize=14, weight='bold')
    ax1.legend(); ax1.grid(True, linestyle='--', alpha=0.3)
    
    ax2.plot(R_vals, implied_ups, 'm-o', label=r'Local $\Upsilon$ Inversion')
    ax2.axhline(0.52, color='k', linestyle='--', label=r'Global Disk Baseline ($\Upsilon=0.52$)')
    ax2.set_xlabel(r'Radius (kpc)'); ax2.set_ylabel(r'Local $\Upsilon$ (Mass/Light)')
    ax2.set_ylim(0, max(implied_ups) + 0.5)
    ax2.legend(); ax2.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, f'IPP_Audit_{galaxy_name}.png')
    plt.savefig(out_path, dpi=300)
    print(f"Audit Saved -> {out_path}")

# ---------------------------------------------------------
# EXECUTION
# ---------------------------------------------------------
if __name__ == "__main__":
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
        
    if os.path.exists(DATA_PATH):
        df = pd.read_csv(DATA_PATH)
        plot_master_correlation(df)
        plot_energy_density_audit(df, 'F571-8')
        plot_energy_density_audit(df, 'D631-7')
    else:
        print(f"Error: {DATA_PATH} not found. Ensure test_sparcV5.py has run.")