import numpy as np
import pandas as pd
import os

original_file = r"C:\zall-pix-iron.fits"
output_csv = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data\desi_168_quasars.csv"

desi_dtype = np.dtype([
    ('TARGETID', '>i8'), ('SURVEY', 'S7'), ('PROGRAM', 'S6'), ('HEALPIX', '>i4'),
    ('SPGRPVAL', '>i4'), ('Z', '>f8'), ('ZERR', '>f8'), ('ZWARN', '>i8'),
    ('CHI2', '>f8'), ('COEFF', '>f8', (10,)), ('NPIXELS', '>i8'),
    ('SPECTYPE', 'S6'), 
    ('SUBTYPE', 'S20'), ('NCOEFF', '>i8'), ('DELTACHI2', '>f8'), 
    ('COADD_FIBERSTATUS', '>i4'), ('TARGET_RA', '>f8'), ('TARGET_DEC', '>f8')
])

ra_min, ra_max = 153.0, 183.0
dec_min, dec_max = -22.0, 8.0

print(f"🔦 Searching for the Deep Laser (Quasars) in the 168° Corridor...")

try:
    file_size = os.path.getsize(original_file)
    header_skip = 250560 
    chunk_size = 800 * 1024 * 1024 
    qso_results = []

    with open(original_file, "rb") as f:
        f.seek(header_skip)
        while f.tell() < file_size:
            raw_bytes = f.read(chunk_size)
            if not raw_bytes: break
            num_rows = len(raw_bytes) // desi_dtype.itemsize
            data = np.frombuffer(raw_bytes, dtype=desi_dtype, count=num_rows)
            
            # BROAD SEARCH: Any Z, but must be SPECTYPE 'QSO'
            stypes = np.array([s.decode('utf-8', errors='ignore').strip('\x00').strip() for s in data['SPECTYPE']])
            
            mask = (data['TARGET_RA'] > ra_min) & (data['TARGET_RA'] < ra_max) & \
                   (data['TARGET_DEC'] > dec_min) & (data['TARGET_DEC'] < dec_max) & \
                   (stypes == 'QSO') & (data['ZWARN'] == 0)
            
            if np.any(mask):
                match = data[mask]
                qso_results.append(pd.DataFrame({
                    'RA': match['TARGET_RA'],
                    'DEC': match['TARGET_DEC'],
                    'Z': match['Z']
                }))
            
            print(f"📉 Search Progress: {(f.tell() / file_size) * 100:.1f}%", end='\r')

    if qso_results:
        final_df = pd.concat(qso_results)
        final_df.to_csv(output_csv, index=False)
        print(f"\n✅ Found {len(final_df)} Quasars in the Corridor!")
    else:
        print("\n❌ Zero Quasars found. This would suggest the 168° axis is even more anomalous than we thought.")

except Exception as e:
    print(f"\n❌ Search failed: {e}")