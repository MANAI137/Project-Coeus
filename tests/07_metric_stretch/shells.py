import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# 1. Path to your Locked Shells
data_dir = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data"

SHELLS = [
    {"file": "shell_1_2mrs_locked_audit.csv", "title": "SHELL 1: The Launch (2MRS)", "color": "indigo"},
    {"file": "shell_2_6dfgs_locked_audit.csv", "title": "SHELL 2: The Boundary (6dFGS)", "color": "darkorange"},
    {"file": "shell_3_sdss_locked_audit.csv",  "title": "SHELL 3: The Bridge (SDSS)", "color": "darkgreen"},
    {"file": "shell_4_desi_locked_audit.csv",  "title": "SHELL 4: The Basin (DESI)", "color": "darkred"}
]

fig = plt.figure(figsize=(20, 18))
bins = np.arange(0, 361, 45) + 15
angles = np.deg2rad(bins[:-1])
width = np.deg2rad(45)

for i, shell in enumerate(SHELLS):
    file_path = os.path.join(data_dir, shell["file"])
    if not os.path.exists(file_path):
        print(f"⚠️ Skipping {shell['file']} (Not found)")
        continue
        
    df = pd.read_csv(file_path)
    df.columns = [c.upper() for c in df.columns]
    total_n = len(df)
    
    # Calculate Sector Counts
    sector_counts = []
    for j in range(len(bins)-1):
        count = len(df[(df['RA'] >= bins[j]) & (df['RA'] < bins[j+1])])
        sector_counts.append(count)
    
    # Create Subplot
    ax = fig.add_subplot(2, 2, i+1, projection='polar')
    bars = ax.bar(angles, sector_counts, width=width, color=shell["color"], alpha=0.6, edgecolor='black')
    
    # Data Overlay
    for j, (angle, count) in enumerate(zip(angles, sector_counts)):
        percent = (count / total_n) * 100
        label = f"N={count}\n({percent:.1f}%)"
        ax.text(angle, count + (max(sector_counts)*0.08), label, ha='center', va='center', fontweight='bold')
        
        # Highlight Metric Core (150-240)
        if 150 <= np.degrees(angle) <= 240:
            bars[j].set_alpha(0.9)
            bars[j].set_linewidth(2)

    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_xticklabels(['15°', '60°', '105°', '150°', '195°', '240°', '285°', '330°'])
    ax.set_title(f"{shell['title']}\nTotal N = {total_n}", fontsize=14, pad=30)

plt.suptitle("PROJECT COEUS: 4-SHELL STRATIGRAPHIC AUDIT\nThe 168° Nozzle to 240° Basin Migration", fontsize=22, y=1.02)
plt.tight_layout()
plt.show()