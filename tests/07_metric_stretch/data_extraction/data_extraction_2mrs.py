import pandas as pd
import numpy as np
import os

# 1. Setup Local Pathing
current_dir = os.path.dirname(os.path.abspath(__file__))
local_data_dir = os.path.join(current_dir, 'data')

input_file = os.path.join(local_data_dir, '2mrs_table3.tsv')
output_csv = os.path.join(local_data_dir, '2mrs_360_bulk_audit.csv')

print(f"🚀 EXECUTING 2MRS 360° BULK AUDIT (Infrared Witness)...")

try:
    # 2. Load and Initial Filter
    # Skipping 68 lines for VizieR standard formatting
    df = pd.read_csv(input_file, sep='|', skiprows=68, header=None, engine='python', on_bad_lines='skip')
    
    # VizieR indices for 2MRS: 4=RA, 5=Dec, 27=cz (velocity), 7=GLAT(b)
    ra_vals = pd.to_numeric(df[4], errors='coerce')
    dec_vals = pd.to_numeric(df[5], errors='coerce')
    cz_vals = pd.to_numeric(df[27], errors='coerce')
    glat_vals = pd.to_numeric(df[7], errors='coerce')

    df_clean = pd.DataFrame({
        'RA': ra_vals, 
        'DEC': dec_vals, 
        'Z': cz_vals / 299792.458, # Convert velocity to Redshift
        'GLAT': glat_vals
    }).dropna()

    # 3. Apply "Bloom Zone" Filters
    # 2MRS is local (high infrared sensitivity). 
    # We audit z=0.01 to 0.1 to capture the "Local Metric Infrastructure".
    z_min, z_max = 0.01, 0.1
    
    # We keep the Galactic Plane Guard (|b| > 10) to avoid noisy data in the disk
    mask = (df_clean['Z'] >= z_min) & \
           (df_clean['Z'] <= z_max) & \
           (np.abs(df_clean['GLAT']) > 10.0)

    final_df = df_clean[mask][['RA', 'DEC', 'Z']]

    # 4. Save the Bulk Population
    if not final_df.empty:
        final_df.to_csv(output_csv, index=False)
        print(f"✅ SUCCESS: {len(final_df)} objects saved for the 360° 2MRS Audit.")
        print(f"📊 Redshift Range: {z_min} - {z_max} (Local Shell)")
        print(f"📁 File: {output_csv}")
    else:
        print("⚠️ No objects matched the criteria. Check VizieR indices.")

except Exception as e:
    print(f"❌ EXTRACTION FAILED: {e}")