import pandas as pd
import numpy as np
import os

# 1. Setup Local Pathing
current_dir = os.path.dirname(os.path.abspath(__file__))
local_data_dir = os.path.join(current_dir, 'data')

input_file = os.path.join(local_data_dir, '6dfgs.tsv')
output_csv = os.path.join(local_data_dir, '6dfgs_360_bulk_audit.csv')

print(f"🔭 EXECUTING 6dFGS 360° BULK AUDIT...")

try:
    # 2. Load Raw TSV (Skipping metadata headers)
    # 6dFGS often has a long header (usually around 50-60 lines)
    df = pd.read_csv(input_file, sep='|', skiprows=56, header=None, engine='python', on_bad_lines='skip')
    
    # Standard 6dFGS column mapping:
    # col 0: RA, col 1: DEC, col 10: cz (velocity), col 13: quality
    ra_vals = pd.to_numeric(df[0], errors='coerce')
    dec_vals = pd.to_numeric(df[1], errors='coerce')
    cz_vals = pd.to_numeric(df[10], errors='coerce')
    quality = pd.to_numeric(df[13], errors='coerce')

    df_clean = pd.DataFrame({
        'RA': ra_vals, 
        'DEC': dec_vals, 
        'Z': cz_vals / 299792.458, # Convert velocity to Redshift
        'Q': quality
    }).dropna()
    
    # 3. Apply Quality & "Bloom Zone" Filters
    # Quality >= 3 (reliable redshifts)
    # Z range: 0.05 to 0.25 (standard 6dF depth)
    z_min, z_max = 0.05, 0.25
    
    mask = (df_clean['Q'] >= 3) & \
           (df_clean['Z'] >= z_min) & \
           (df_clean['Z'] <= z_max) & \
           (df_clean['DEC'] < 0) # Strictly Southern Hemisphere

    final_df = df_clean[mask][['RA', 'DEC', 'Z']]

    # 4. Save the Bulk Population
    if not final_df.empty:
        final_df.to_csv(output_csv, index=False)
        print(f"✅ SUCCESS: {len(final_df)} objects saved for the 360° Southern Audit.")
        print(f"📊 Redshift Range: {z_min} - {z_max}")
        print(f"📁 File: {output_csv}")
    else:
        print("⚠️ No objects matched the criteria. Check your .tsv header skip count.")

except Exception as e:
    print(f"❌ EXTRACTION FAILED: {e}")