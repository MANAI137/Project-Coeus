import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load our two populations
df_gal = pd.read_csv(r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data\desi_168_spectral.csv")
df_qso = pd.read_csv(r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data\desi_168_quasars.csv")

# Filter for the Galaxy population (mostly low Z)
df_gal = df_gal[df_gal['SPECTYPE'].str.strip() == 'GALAXY']

def get_stats(data, bins):
    w, e, z_c = [], [], []
    for i in range(len(bins)-1):
        shell = data[(data['Z'] >= bins[i]) & (data['Z'] < bins[i+1])]
        if len(shell) > 10:
            spread = np.std(shell['RA'])
            err = spread / np.sqrt(2 * (len(shell) - 1))
            w.append(spread)
            e.append(err)
            z_c.append((bins[i] + bins[i+1]) / 2)
    return z_c, w, e

# Broad Z-range to capture both the "Bloom" and the "Deep Laser"
z_bins = np.arange(0.05, 3.5, 0.2)

z_gal, w_gal, e_gal = get_stats(df_gal, z_bins)
z_qso, w_qso, e_qso = get_stats(df_qso, z_bins)

plt.figure(figsize=(12, 7))

# The Bloom (Galaxies)
plt.errorbar(z_gal, w_gal, yerr=e_gal, fmt='s-', color='orange', 
             label='Galaxies (Low Energy / The Bloom)', capsize=4)

# The Laser (Quasars)
plt.errorbar(z_qso, w_qso, yerr=e_qso, fmt='o-', color='purple', 
             label='Quasars (High Energy / The Laser)', capsize=4, linewidth=2)

plt.title("The Refractive Audit: Do Quasars stay in the 'Bore'?")
plt.xlabel("Redshift (z)")
plt.ylabel("Angular Aperture (RA Std Dev)")
plt.gca().invert_xaxis()
plt.legend()
plt.grid(True, alpha=0.2)
plt.show()