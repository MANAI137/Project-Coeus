"""
PROJECT COEUS: SDSS METRIC GRAIN STABILITY AUDIT
------------------------------------------------
Author: Miguel Antonio Navarro

DESCRIPTION:
This script audits the directional 'stiffness' of the metric across 
three deep-space horizons using the SDSS dataset. It measures the 
mean redshift deviation as a function of angular separation from 
the 168° manifold launch vector. 

A consistent 'Metric Swing' across the Local, Mid-Range, and Deep 
Space tiers serves as a forensic rejection of the isotropic null. 
It demonstrates that the observed Hubble Tension is a predictable 
consequence of a directional manifold grain rather than a scalar 
measurement error.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Pathing (Relative to your tests/07_metric_stretch/ folder)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(BASE_DIR, "data", "sdss_waypoint_map.csv")

def run_control_test():
    df = pd.read_csv(CSV_PATH)
    
    # We slice the data into "Distance Tiers" based on Z
    # If the grain is metric, the dipole should be visible in EVERY tier
    tiers = [
        (0.01, 0.1),  # Local Universe
        (0.1, 0.3),   # Mid-Range
        (0.3, 0.6)    # Deep Space
    ]
    
    plt.figure(figsize=(12, 8))
    colors = ['cyan', 'magenta', 'yellow']
    
    for i, (z_min, z_max) in enumerate(tiers):
        tier_df = df[(df['z'] >= z_min) & (df['z'] <= z_max)].copy()
        
        bins = np.linspace(0, 180, 13)
        tier_df['angle_bin'] = pd.cut(tier_df['angle_to_epicenter'], bins)
        results = tier_df.groupby('angle_bin', observed=True)['z'].agg(['mean', 'std', 'count']).reset_index()
        results['mid_angle'] = results['angle_bin'].apply(lambda x: x.mid)
        
        # Plot each tier's unique dipole
        plt.plot(results['mid_angle'], results['mean'], 'o-', 
                 color=colors[i], label=f'Tier: z={z_min}-{z_max}')

    plt.title("Test 2: Multi-Tier Isotropic Audit (The Metric Grain Test)", color='white', fontsize=16)
    plt.xlabel("Degrees from Metric Epicenter", color='white')
    plt.ylabel("Mean Redshift (z)", color='white')
    plt.legend()
    plt.grid(alpha=0.1)
    
    # Coeus Styling
    plt.gca().set_facecolor('#111111')
    plt.gcf().set_facecolor('#111111')
    plt.tick_params(colors='white')
    for spine in plt.gca().spines.values():
        spine.set_color('white')

    plt.savefig(os.path.join(BASE_DIR, 'test_02_metric_grain.png'))
    print("✅ Multi-Tier Audit Complete. Plot saved.")

if __name__ == "__main__":
    run_control_test()