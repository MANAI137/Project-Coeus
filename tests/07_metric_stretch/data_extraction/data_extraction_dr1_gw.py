import numpy as np
import pandas as pd
import os

original_file = r"C:\zall-pix-iron.fits"
# Renamed output to reflect the Basin Audit
output_csv = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\desi_240_gw_audit.csv"

# The validated DESI Iron Row Structure remains the same
desi_dtype = np.dtype([
    ('TARGETID', '>i8'), ('SURVEY', 'S7'), ('PROGRAM', 'S6'), ('HEALPIX', '>i4'),
    ('SPGRPVAL', '>i4'), ('Z', '>f8'), ('ZERR', '>f8'), ('ZWARN', '>i8'),
    ('CHI2', '>f8'), ('COEFF', '>f8', (10,)), ('NPIXELS', '>i8'),
    ('SPECTYPE', 'S6'), ('SUBTYPE', 'S20'), ('NCOEFF', '>i8'),
    ('DELTACHI2', '>f8'), ('COADD_FIBERSTATUS', '>i4'),
    ('TARGET_RA', '>f8'), ('TARGET_DEC', '>f8')
])

# --- TARGET: THE 240° BASIN AUDIT ---
# Centered on 240, +/- 30 degrees RA to see the "Delta"
ra_min, ra_max = 210.0, 270.0 
# Expanded DEC to +/- 40 to test North/South Time Symmetry
dec_min, dec_max = -40.0, 40.0

print(f"🚀 Starting Basin Audit Streamer on {original_file}...")
print(f"Targeting RA [{ra_min}-{ra_max}] | DEC [{dec_min}-{dec_max}]")

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
            if num_rows == 0: break
            
            # Parse the chunk
            data = np.frombuffer(raw_bytes, dtype=desi_dtype, count=num_rows)
            
            # Filter Logic: Targeting the Basin + Quality Control
            mask = (data['TARGET_RA'] >= ra_min) & (data['TARGET_RA'] <= ra_max) & \
                   (data['TARGET_DEC'] >= dec_min) & (data['TARGET_DEC'] <= dec_max) & \
                   (data['ZWARN'] == 0) & \
                   (data['Z'] > 0.001) & (data['Z'] < 0.6) # Focus on structural maturity range
            
            if np.any(mask):
                match = data[mask]
                all_results.append(pd.DataFrame({
                    'RA': match['TARGET_RA'],
                    'DEC': match['TARGET_DEC'],
                    'Z': match['Z']
                }))
            
            progress = (f.tell() / file_size) * 100
            print(f"📉 Progress: {progress:.1f}% | Basin Count: {sum(len(x) for x in all_results)}", end='\r')

    # Save the full Basin dataset
    if all_results:
        final_df = pd.concat(all_results)
        final_df.to_csv(output_csv, index=False)
        print(f"\n✅ Success! Saved {len(final_df)} objects to {output_csv}")
    else:
        print("\n❌ No objects found in the target Basin coordinates.")

except Exception as e:
    print(f"\n❌ Streaming failed: {e}")