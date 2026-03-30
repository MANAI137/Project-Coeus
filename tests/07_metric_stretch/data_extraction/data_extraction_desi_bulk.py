import numpy as np
import pandas as pd
import os

# Paths
original_file = r"C:\zall-pix-iron.fits"
output_csv = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data\desi_360_bulk_audit.csv"

desi_dtype = np.dtype([
    ('TARGETID', '>i8'), ('SURVEY', 'S7'), ('PROGRAM', 'S6'), ('HEALPIX', '>i4'),
    ('SPGRPVAL', '>i4'), ('Z', '>f8'), ('ZERR', '>f8'), ('ZWARN', '>i8'),
    ('CHI2', '>f8'), ('COEFF', '>f8', (10,)), ('NPIXELS', '>i8'),
    ('SPECTYPE', 'S6'), 
    ('SUBTYPE', 'S20'), ('NCOEFF', '>i8'), ('DELTACHI2', '>f8'), 
    ('COADD_FIBERSTATUS', '>i4'), ('TARGET_RA', '>f8'), ('TARGET_DEC', '>f8')
])

# FULL CIRCLE SCAN: RA 0 to 360
# Keeping the High-DEC "Clear Window" (10° to 40°)
ra_min, ra_max = 0.0, 360.0
dec_min, dec_max = 10.0, 40.0 

print(f"🔭 Executing Full 360° MFA Scan (High-DEC: 10° to 40°)...")

try:
    file_size = os.path.getsize(original_file)
    header_skip = 250560 
    chunk_size = 1024 * 1024 * 1024 # 1GB chunks for faster throughput
    results = []

    with open(original_file, "rb") as f:
        f.seek(header_skip)
        while f.tell() < file_size:
            raw_bytes = f.read(chunk_size)
            if not raw_bytes: break
            num_rows = len(raw_bytes) // desi_dtype.itemsize
            data = np.frombuffer(raw_bytes, dtype=desi_dtype, count=num_rows)
            
            # Filtering for quality, redshift, and the 360° High-DEC window
            mask = (data['TARGET_RA'] >= ra_min) & (data['TARGET_RA'] <= ra_max) & \
                   (data['TARGET_DEC'] >= dec_min) & (data['TARGET_DEC'] <= dec_max) & \
                   (data['ZWARN'] == 0) & (data['Z'] > 0.001)
            
            if np.any(mask):
                match = data[mask]
                results.append(pd.DataFrame({
                    'RA': match['TARGET_RA'], 
                    'DEC': match['TARGET_DEC'], 
                    'Z': match['Z']
                }))
            
            print(f"📉 360° Scan Progress: {(f.tell() / file_size) * 100:.1f}%", end='\r')

    if results:
        final_df = pd.concat(results)
        final_df.to_csv(output_csv, index=False)
        print(f"\n✅ All-Sky MFA Scan Complete! Found {len(final_df)} objects.")
        print(f"📁 Data saved to: {output_csv}")
    else:
        print("\n❌ No data found in the 360° High-DEC window.")

except Exception as e:
    print(f"\n❌ 360° Scan failed: {e}")