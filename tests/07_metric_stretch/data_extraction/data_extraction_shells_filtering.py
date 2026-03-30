import pandas as pd
import os

# 1. Update these to your actual bulk file names
data_dir = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data"

LOCKDOWN_CONFIG = {
    "2mrs_360_bulk_audit.csv": {"z_min": 0.010, "z_max": 0.045, "name": "SHELL_1_2MRS"},
    "6dfgs_360_bulk_audit.csv": {"z_min": 0.060, "z_max": 0.110, "name": "SHELL_2_6DFGS"},
    "sdss_360_bulk_audit.csv":  {"z_min": 0.130, "z_max": 0.180, "name": "SHELL_3_SDSS"},
    "desi_360_bulk_audit.csv":  {"z_min": 0.220, "z_max": 0.350, "name": "SHELL_4_DESI"}
}

def safe_filter():
    print("🛡️ STARTING MEMORY-SAFE STRATIGRAPHIC FILTER...")
    for filename, limits in LOCKDOWN_CONFIG.items():
        file_path = os.path.join(data_dir, filename)
        output_name = os.path.join(data_dir, f"{limits['name']}_locked_audit.csv")
        
        if not os.path.exists(file_path):
            print(f"⚠️ Skipping missing file: {filename}")
            continue

        print(f"Processing {filename}...")
        
        # Open output file and write header
        first_chunk = True
        total_count = 0
        
        # CHUNK PROCESSING (The Crash-Proof Part)
        for chunk in pd.read_csv(file_path, chunksize=50000):
            chunk.columns = [c.upper() for c in chunk.columns]
            
            # Apply the Z-Gap mask
            mask = (chunk['Z'] >= limits['z_min']) & (chunk['Z'] <= limits['z_max'])
            filtered_chunk = chunk[mask]
            
            # Write to disk immediately
            mode = 'w' if first_chunk else 'a'
            header = True if first_chunk else False
            filtered_chunk.to_csv(output_name, mode=mode, header=header, index=False)
            
            first_chunk = False
            total_count += len(filtered_chunk)
            
        print(f"✅ Created {limits['name']} | Final N={total_count}")

if __name__ == "__main__":
    safe_filter()