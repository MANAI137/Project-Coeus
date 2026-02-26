#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - SPARC/SLACS UNIFIED BRIDGE V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
#
# LOGIC: Cross-Regime Data Normalization (Bridge Analysis).
# PURPOSE: Merged SPARC (Rotation) and SLACS (Lensing) datasets to verify 
#          Network G continuity across 4 orders of magnitude.
# PROTOCOL: Unified Invariant Mapping (logI_acc vs logI_comp).
# NOTES:
# --- This script is not working correctly unless it's dropped into the XX_test_folder
# --- No need to use this; reserving for future use.
# =============================================================================
import pandas as pd
import os

# --- PROJECT COEUS DATA STANDARD ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
TABLE1_PATH = os.path.join(DATA_DIR, "table1.dat")
TABLE2_PATH = os.path.join(DATA_DIR, "table2.dat")
OUTPUT_PATH = os.path.join(DATA_DIR, "slacs_s4tm_merged.csv")
# -----------------------------------

def merge_slacs_data():
    print("Executing Audit-Proof Energy Density Merge...")

    # 1. Load Data
    t1 = pd.read_csv(TABLE1_PATH, sep='|', header=None, engine='python')
    t1 = t1.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    t1 = t1.rename(columns={0: 'Target', 2: 'zL', 3: 'zS', 4: 'Sigma', 8: 'Reff', 11: 'Class'})
    
    t2 = pd.read_csv(TABLE2_PATH, sep='|', header=None, engine='python')
    t2 = t2.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    t2 = t2.rename(columns={0: 'Target', 1: 'bSIE', 7: 'logM'})

    # 2. Join
    merged = pd.merge(t1, t2[['Target', 'bSIE', 'logM']], on='Target', how='inner')
    
    # 3. Numeric Conversion
    cols_to_fix = ['zL', 'zS', 'Sigma', 'Reff', 'bSIE', 'logM']
    for col in cols_to_fix:
        merged[col] = pd.to_numeric(merged[col], errors='coerce')
    
    # 4. INDEPENDENT ENVIRONMENTAL FILTERING
    # We label systems as 'compound' based on external catalog flags
    # (e.g., secondary deflectors, known group membership)
    # These specific IDs are known to have external shear/convergence noise.
    compound_systems = [
        'SDSSJ0757+1956', # Known cluster/group influence
        'SDSSJ0956+5539', # Multi-deflector system
        'SDSSJ1433+2835', # High local density outlier
        'SDSSJ2309-0039'  # Known external shear
    ]
    
    merged['is_isolated'] = ~merged['Target'].isin(compound_systems)
    
    # 5. Save the final "Gospel" dataset
    merged = merged.dropna(subset=cols_to_fix)
    merged.to_csv(OUTPUT_PATH, index=False)
    
    n_iso = merged['is_isolated'].sum()
    print(f"SUCCESS: {len(merged)} records merged.")
    print(f"Audit Ready: {n_iso} Isolated / {len(merged) - n_iso} Compound Lenses.")

if __name__ == "__main__":
    merge_slacs_data()