import pandas as pd
import numpy as np

# Load data
qso_path = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data\dr16q_lean_audit.csv"
df = pd.read_csv(qso_path)

# Theoretical Vectors
STRIKE_ORIGIN = 168.0
BASIN_CENTER = 239.1

def calculate_alignment_score(ra_array, target_ra):
    # Measures how "bunched" the RA is around a specific target
    # 0 = Random, 1 = Perfectly Aligned
    diff = np.abs(ra_array - target_ra)
    diff = np.where(diff > 180, 360 - diff, diff)
    return 1 - (np.mean(diff) / 90.0)

# Audit across Shells
z_bins = [0.1, 1.0, 2.0, 3.0, 5.0]
print("--- PROJECT COEUS: CMB ALIGNMENT AUDIT ---")
for i in range(len(z_bins)-1):
    z_min, z_max = z_bins[i], z_bins[i+1]
    shell = df[(df['Z'] >= z_min) & (df['Z'] <= z_max)]
    
    strike_score = calculate_alignment_score(shell['RA'], STRIKE_ORIGIN)
    basin_score = calculate_alignment_score(shell['RA'], BASIN_CENTER)
    
    print(f"Shell Z[{z_min}-{z_max}]:")
    print(f"  Alignment with Strike (168°): {strike_score:.4f}")
    print(f"  Alignment with Basin  (240°): {basin_score:.4f}")