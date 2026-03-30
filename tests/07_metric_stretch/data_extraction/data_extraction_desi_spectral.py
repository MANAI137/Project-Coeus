import numpy as np
import pandas as pd
import os

original_file = r"C:\zall-pix-iron.fits"
output_csv = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\desi_168_spectral.csv"

# The DESI Iron Row Structure with SPECTYPE included
desi_dtype = np.dtype([
    ('TARGETID', '>i8'), ('SURVEY', 'S7'), ('PROGRAM', 'S6'), ('HEALPIX', '>i4'),
    ('SPGRPVAL', '>i4'), ('Z', '>f8'), ('ZERR', '>f8'), ('ZWARN', '>i8'),
    ('CHI2', '>f8'), ('COEFF', '>f8', (10,)), ('NPIXELS', '>i8'),
    ('SPECTYPE', 'S6'), # <--- This is the key we were missing
    ('SUBTYPE', 'S20'), ('NCOEFF', '>i8'), ('DELTACHI2', '>f8'), 
    ('COADD_FIBERSTATUS', '>i4'), ('TARGET_RA', '>f8'), ('TARGET_DEC', '>f8')
])

ra_min, ra_max = 153.0, 183.0
dec_min, dec_max = -22.0, 8.0

print(f"🚀 Re-streaming 168° Corridor with Spectral Classification...")

try:
    file_size = os.path.getsize(original_file)
    header_skip = 250560 
    chunk_size = 500 * 1024 * 1024 
    all_results = []

    with open(original_file, "rb") as f:
        f.seek(header_skip)
        while f.tell() < file_size:
            raw_bytes = f.read(chunk_size)
            if not raw_bytes: break
            num_rows = len(raw_bytes) // desi_dtype.itemsize
            data = np.frombuffer(raw_bytes, dtype=desi_dtype, count=num_rows)
            
            mask = (data['TARGET_RA'] > ra_min) & (data['TARGET_RA'] < ra_max) & \
                   (data['TARGET_DEC'] > dec_min) & (data['TARGET_DEC'] < dec_max) & \
                   (data['ZWARN'] == 0) & (data['Z'] > 0.001)
            
            if np.any(mask):
                match = data[mask]
                all_results.append(pd.DataFrame({
                    'RA': match['TARGET_RA'],
                    'DEC': match['TARGET_DEC'],
                    'Z': match['Z'],
                    'SPECTYPE': match['SPECTYPE'].astype(str).str.strip() # Clean the string
                }))
            print(f"📉 Progress: {(f.tell() / file_size) * 100:.1f}%", end='\r')

    final_df = pd.concat(all_results)
    final_df.to_csv(output_csv, index=False)
    print(f"\n✅ Created spectral dataset: {len(final_df)} entries.")

except Exception as e:
    print(f"\n❌ Streaming failed: {e}")