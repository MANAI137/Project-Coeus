import numpy as np
import pandas as pd
import os

original_file = r"C:\zall-pix-iron.fits"
output_csv = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\desi_168_corridor.csv"

# The validated DESI Iron Row Structure
desi_dtype = np.dtype([
    ('TARGETID', '>i8'), ('SURVEY', 'S7'), ('PROGRAM', 'S6'), ('HEALPIX', '>i4'),
    ('SPGRPVAL', '>i4'), ('Z', '>f8'), ('ZERR', '>f8'), ('ZWARN', '>i8'),
    ('CHI2', '>f8'), ('COEFF', '>f8', (10,)), ('NPIXELS', '>i8'),
    ('SPECTYPE', 'S6'), ('SUBTYPE', 'S20'), ('NCOEFF', '>i8'),
    ('DELTACHI2', '>f8'), ('COADD_FIBERSTATUS', '>i4'),
    ('TARGET_RA', '>f8'), ('TARGET_DEC', '>f8')
])

# 168° Corridor Boundaries
ra_min, ra_max = 153.0, 183.0
dec_min, dec_max = -22.0, 8.0

print(f"🚀 Starting Full Streamer on {original_file}...")

try:
    file_size = os.path.getsize(original_file)
    header_skip = 250560 # The offset we validated
    chunk_size = 500 * 1024 * 1024 # 500MB chunks for speed/stability
    
    all_results = []

    with open(original_file, "rb") as f:
        f.seek(header_skip)
        
        while f.tell() < file_size:
            raw_bytes = f.read(chunk_size)
            if not raw_bytes: break
            
            num_rows = len(raw_bytes) // desi_dtype.itemsize
            if num_rows == 0: break
            
            # Parse the chunk
            data = np.frombuffer(raw_bytes, dtype=desi_dtype, count=num_rows)
            
            # Filter
            mask = (data['TARGET_RA'] > ra_min) & (data['TARGET_RA'] < ra_max) & \
                   (data['TARGET_DEC'] > -22.0) & (data['TARGET_DEC'] < 8.0) & \
                   (data['ZWARN'] == 0) & (data['Z'] > 0.001) # Filter out the z=0 noise
            
            if np.any(mask):
                match = data[mask]
                all_results.append(pd.DataFrame({
                    'RA': match['TARGET_RA'],
                    'DEC': match['TARGET_DEC'],
                    'Z': match['Z']
                }))
            
            progress = (f.tell() / file_size) * 100
            print(f"📉 Progress: {progress:.1f}% | Corridor Count: {sum(len(x) for x in all_results)}", end='\r')

    # Save the full corridor dataset
    if all_results:
        final_df = pd.concat(all_results)
        final_df.to_csv(output_csv, index=False)
        print(f"\n✅ Success! Saved {len(final_df)} galaxies to {output_csv}")
    else:
        print("\n❌ No galaxies found beyond the z=0 threshold.")

except Exception as e:
    print(f"\n❌ Streaming failed: {e}")