import astropy.io.fits as fits
import pandas as pd
import numpy as np

# 1. Configuration
file_path = r'D:\specObj-dr17.fits'
output_csv = 'sdss_360_bulk_audit.csv'

# Metric Epicenter (CMB Dipole Reference)
RA_DIPOLE = 168.0
DEC_DIPOLE = -7.0

def angular_separation(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    a = np.sin((dec2-dec1)/2)**2 + np.cos(dec1)*np.cos(dec2)*np.sin((ra2-ra1)/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))

print(f"🚀 Processing SDSS DR17 Bulk Flow...")

with fits.open(file_path, memmap=True) as hdul:
    data = hdul[1].data
    
    # 2. String-Resistant Masking
    # We strip spaces and ensure case-insensitivity for 'GALAXY'
    # We also filter for high-confidence redshifts (ZWARNING == 0)
    print("🔭 Filtering 5.8M rows (Stripping strings and checking quality)...")
    
    # Converting to string to handle byte-literals and padding
    classes = np.char.strip(data['CLASS'].astype(str))
    mask = (classes == 'GALAXY') & (data['ZWARNING'] == 0) & (data['Z'] > 0.01)
    
    if not np.any(mask):
        print("❌ Mask failed. Checking first few entries of 'CLASS' manually...")
        print(f"Sample data['CLASS'][:5]: {data['CLASS'][:5]}")
    else:
        ra = data['PLUG_RA'][mask]
        dec = data['PLUG_DEC'][mask]
        z = data['Z'][mask]
        
        print(f"✅ Found {len(z)} Galaxies. Calculating Vectors...")
        dist_from_epicenter = angular_separation(ra, dec, RA_DIPOLE, DEC_DIPOLE)

        # 3. Create Bulk Audit DataFrame
        df = pd.DataFrame({
            'RA': ra,
            'DEC': dec,
            'Z': z,
            'angle_to_epicenter': dist_from_epicenter
        })

        # Save the full result
        print(f"💾 Saving to {output_csv}...")
        df.to_csv(output_csv, index=False)
        print("✅ SDSS Bulk Flow extraction complete.")