"""
PROJECT COEUS: MANIFOLD BORE & EXPANSION AUDIT
----------------------------------------------
Author: Miguel Antonio Navarro
Version: 1.0 (March 2026)

DESCRIPTION:
This script performs a deep-time geometric analysis of the 168°–240° manifold. 
By calculating the Angular Standard Deviation (the 'Bore') of Quasar populations 
across increasing redshift shells (z=0.5 to 4.5), it tests for the existence 
of a collimated structural corridor.

Unlike a simple density count, this audit measures the geometric stability 
of the metric itself. It compares the 240° Basin (the downstream target) 
against its 60° antipode (the isotropic control). A persistent, narrower 
spread in the Basin vs. the Control functions as a 'Metric Fingerprint,' 
proving that the observed anisotropy is a stable, long-lived feature 
of the cosmic vacuum rather than a local density fluctuation.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# 1. Path Configuration
data_dir = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data"
qso_path = os.path.join(data_dir, "dr16q_lean_audit.csv")

def run_geometry_audit():
    if not os.path.exists(qso_path):
        print(f"❌ File not found: {qso_path}")
        return

    # Load Data
    df = pd.read_csv(qso_path)
    df.columns = [c.upper() for c in df.columns]

    # 2. Define the Sectors
    # Basin (Target) vs Control (Opposite Pole)
    sectors = {
        "Basin (240°)": {"ra_min": 210, "ra_max": 270, "color": "red"},
        "Control (60°)": {"ra_min": 30, "ra_max": 90, "color": "blue"}
    }

    # 3. Define Redshift Shells
    z_bins = np.arange(0.5, 4.5, 0.5)
    results = []

    print(f"📐 Starting Geometric Spread Audit across {len(z_bins)-1} shells...")

    for i in range(len(z_bins)-1):
        z_min, z_max = z_bins[i], z_bins[i+1]
        z_mid = (z_min + z_max) / 2
        
        shell_df = df[(df['Z'] >= z_min) & (df['Z'] <= z_max)]
        
        for name, bounds in sectors.items():
            # Extract the sector
            mask = (shell_df['RA'] >= bounds['ra_min']) & (shell_df['RA'] <= bounds['ra_max'])
            sector_data = shell_df[mask]
            
            if len(sector_data) > 50:
                # Calculate Angular Spread (Standard Deviation of RA)
                spread = sector_data['RA'].std()
                results.append({
                    "Z": z_mid,
                    "Sector": name,
                    "Spread": spread,
                    "Count": len(sector_data),
                    "Color": bounds['color']
                })

    res_df = pd.DataFrame(results)

    # 4. Visualization: The "Cone of the Spray"
    plt.figure(figsize=(12, 7))
    
    for name, group in res_df.groupby("Sector"):
        plt.plot(group['Z'], group['Spread'], marker='o', label=name, 
                 color=group['Color'].iloc[0], linewidth=2.5)
        
        # Add Count Labels
        for x, y, n in zip(group['Z'], group['Spread'], group['Count']):
            plt.text(x, y + 0.1, f"n={n}", ha='center', fontsize=9, alpha=0.7)

    plt.axhline(17.3, color='gray', linestyle='--', alpha=0.5, label='Isotropic Expected (Random)')
    
    plt.title("PROJECT COEUS: GEOMETRIC SPREAD AUDIT\nAngular Width (RA Std) vs. Cosmic Depth (Z)", fontsize=16)
    plt.xlabel("Redshift (Z) - Deep Time Horizon →", fontsize=12)
    plt.ylabel("Angular Spread (Degrees σ)", fontsize=12)
    plt.grid(alpha=0.3)
    plt.legend()
    
    # Forensic Indicator
    plt.annotate('Launch Origin memory?', xy=(4.0, 16.5), xytext=(3.0, 15.5),
                 arrowprops=dict(facecolor='black', shrink=0.05))

    plt.show()

    # 5. Summary Report
    print("\n--- GEOMETRIC SPREAD SUMMARY ---")
    pivot = res_df.pivot(index="Z", columns="Sector", values="Spread")
    print(pivot)

if __name__ == "__main__":
    run_geometry_audit()