"""
PROJECT COEUS: METRIC GRAIN STABILITY AUDIT
-------------------------------------------
Author: Miguel Antonio Navarro

DESCRIPTION:
This script tests the directional 'stiffness' of the cosmic manifold. 
It measures the Mean Observed Redshift (z) as a function of angular 
separation from the 168° launch vector. 

By comparing three distinct distance horizons (Local, Mid, and Deep), 
it determines if the 'Metric Swing' (anisotropic deviation) is a 
persistent feature of the vacuum geometry or merely a local cluster 
artifact. A consistent swing across all tiers indicates a global 
directional grain in the metric itself, challenging the isotropic null.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Pathing (Directly targeting your processed SDSS set)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(BASE_DIR, "data", "2mrs_standalone_buoys.csv") # Updated for consistency

def run_multi_tier_audit():
    if not os.path.exists(CSV_PATH):
        print(f"❌ ERROR: Could not find {CSV_PATH}")
        return

    df = pd.read_csv(CSV_PATH)
    print(f"📡 Analyzing {len(df)} waypoints across 3 horizons...")

    # Distance Tiers
    tiers = [
        (0.01, 0.1),  # Local Universe (matches 2MRS scale)
        (0.1, 0.3),   # Mid-Range
        (0.3, 0.6)    # Deep Space (Metric maturity)
    ]
    
    plt.figure(figsize=(12, 8))
    # Using the Project Coeus color palette
    colors = ['#00FFFF', '#FF00FF', '#FFFF00'] # Cyan, Magenta, Yellow
    
    for i, (z_min, z_max) in enumerate(tiers):
        tier_df = df[(df['z'] >= z_min) & (df['z'] <= z_max)].copy()
        
        if tier_df.empty:
            print(f"⚠️ Tier z={z_min}-{z_max} is empty. Skipping.")
            continue

        bins = np.linspace(0, 180, 13)
        tier_df['angle_bin'] = pd.cut(tier_df['angle_to_epicenter'], bins)
        
        results = tier_df.groupby('angle_bin', observed=True)['z'].agg(['mean', 'std', 'count']).reset_index()
        results['mid_angle'] = results['angle_bin'].apply(lambda x: x.mid)
        
        # Calculate the relative swing for the ledger
        swing = ((results['mean'].max() - results['mean'].min()) / results['mean'].mean()) * 100
        print(f"📊 Tier z={z_min}-{z_max}: {len(tier_df)} galaxies | Metric Swing: {swing:.2f}%")

        # Plot each tier
        plt.plot(results['mid_angle'], results['mean'], 'o-', 
                 color=colors[i], linewidth=2.5, markersize=8,
                 label=f'Tier: z={z_min}-{z_max} ({swing:.1f}% swing)')

    plt.title("Test 2: Multi-Tier Isotropic Audit (The Metric Grain Test)", color='white', fontsize=18, pad=20)
    plt.xlabel("Degrees from Metric Epicenter (168°, -7°)", color='white', fontsize=14)
    plt.ylabel("Mean Observed Redshift (z)", color='white', fontsize=14)
    
    # Legend styling
    leg = plt.legend(facecolor='#111111', edgecolor='white', fontsize=12)
    for text in leg.get_texts():
        text.set_color('white')
        
    plt.grid(color='white', alpha=0.1, linestyle='--')
    
    # Coeus Styling
    plt.gca().set_facecolor('#111111')
    plt.gcf().set_facecolor('#111111')
    plt.tick_params(colors='white', which='both', labelsize=12)
    for spine in plt.gca().spines.values():
        spine.set_color('white')

    save_path = os.path.join(BASE_DIR, 'test_02_metric_grain.png')
    plt.savefig(save_path, dpi=300)
    print(f"✅ Multi-Tier Audit Complete. Plot saved: {save_path}")

if __name__ == "__main__":
    run_multi_tier_audit()