import pandas as pd
import os

# 1. Paths
source_path = r"C:\DR16Q_v4.tsv"
output_csv = r"C:\Users\Miguel\Documents\github_repositories\Project-Coeus\tests\07_metric_stretch\data\dr16q_lean_audit.csv"

def extract_dr16q_surgical():
    print(f"🔭 Searching for data start in {source_path}...")
    
    header_line_index = None
    delimiter = '|' 
    
    with open(source_path, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            if 'RAJ2000' in line and not line.startswith('#'):
                header_line_index = i
                break
                
    if header_line_index is None:
        print("❌ Could not find the header row (RAJ2000) in the file.")
        return

    print(f"🎯 Found header at line {header_line_index + 1}")

    try:
        # Load the data
        df = pd.read_csv(
            source_path,
            sep=delimiter,
            skiprows=header_line_index,
            on_bad_lines='skip',
            low_memory=False
        )
        
        # Clean column names and remove those empty 'Unnamed' columns from the pipes
        df.columns = [str(c).strip() for c in df.columns]
        df = df.loc[:, ~df.columns.str.contains('^Unnamed|^$')]
        
        # Mapping exactly as found in your last run
        # Note: 'z' is lowercase in your file!
        ra_col = next((c for c in df.columns if 'RAJ2000' in c.upper()), None)
        dec_col = next((c for c in df.columns if 'DEJ2000' in c.upper()), None)
        z_col = next((c for c in df.columns if c == 'z' or c == 'Z'), None)

        if not ra_col or not dec_col or not z_col:
            print(f"❌ Column mapping failed. Found columns: {list(df.columns[:5])}")
            return

        print(f"🚀 Slicing data via columns: {ra_col}, {dec_col}, {z_col}")
        
        # Create a new lean dataframe with just our 3 target columns
        df_lean = pd.DataFrame()
        df_lean['RA'] = pd.to_numeric(df[ra_col], errors='coerce')
        df_lean['DEC'] = pd.to_numeric(df[dec_col], errors='coerce')
        df_lean['Z'] = pd.to_numeric(df[z_col], errors='coerce')

        # Drop the "Units" and "Dashes" rows that VizieR includes
        initial_count = len(df_lean)
        df_lean = df_lean.dropna()
        
        print(f"💾 Saving {len(df_lean)} Quasars to CSV.")
        df_lean.to_csv(output_csv, index=False)
        print("✨ Done! The 'High-Energy map' is officially ready.")

    except Exception as e:
        print(f"❌ Extraction failed: {e}")

if __name__ == "__main__":
    extract_dr16q_surgical()