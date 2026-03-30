import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the new spectral dataset
df = pd.read_csv(r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\desi_168_spectral.csv")

def get_stats(data, bins, label):
    # Strip any potential whitespace just in case
    subset = data[data['SPECTYPE'].str.strip() == label]
    z_c, w, e = [], [], []
    for i in range(len(bins)-1):
        shell = subset[(subset['Z'] >= bins[i]) & (subset['Z'] < bins[i+1])]
        if len(shell) > 15: # Ensuring statistical significance
            spread = np.std(shell['RA'])
            err = spread / np.sqrt(2 * (len(shell) - 1))
            w.append(spread)
            e.append(err)
            z_c.append((bins[i] + bins[i+1]) / 2)
    return z_c, w, e

z_bins = np.arange(0.05, 0.45, 0.04)

z_qso, w_qso, e_qso = get_stats(df, z_bins, 'QSO')
z_gal, w_gal, e_gal = get_stats(df, z_bins, 'GALAXY')

plt.figure(figsize=(12, 7))

# The "Refractive" Comparison
plt.errorbar(z_qso, w_qso, yerr=e_qso, fmt='o-', color='purple', label='Quasars (Relativistic Core)', capsize=3, linewidth=2)
plt.errorbar(z_gal, w_gal, yerr=e_gal, fmt='s-', color='orange', label='Galaxies (Metric Spray)', capsize=3, linewidth=2)

plt.title("The Refractive Audit: Does the 168° Vector Sort Matter by Energy?")
plt.xlabel("Redshift (z)")
plt.ylabel("Angular Aperture (RA Std Dev)")
plt.gca().invert_xaxis()
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

# Quick print of the sample sizes for the skeptic (ChatGPT)
print(f"📊 Total Quasars in corridor: {len(df[df['SPECTYPE']=='QSO'])}")
print(f"📊 Total Galaxies in corridor: {len(df[df['SPECTYPE']=='GALAXY'])}")