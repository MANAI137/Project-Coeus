#!/usr/bin/env python3
# # =============================================================================
# PROJECT COEUS - THE FINAL AUDIT V.1.1 (STRICT BLIND)
# Author: Miguel Navarro
# Logic: SI-Grounded Phase Threshold (287.35 M_sun/pc^2).
# Protocol: 25% Anchor / 75% Blind Prediction (REPLICATION-GRADE).
# =============================================================================
import numpy as np
import pandas as pd
import io, os, sys

# 1. ENGINE LINKING
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(BASE_DIR, '../../'))
sys.path.insert(0, PROJECT_ROOT)

try:
    from network_g_core import engineV6 as engine
    print(f"✅ SUCCESS: Network G Engine V5 linked. Initiating Strict Blind Audit.")
except ImportError:
    sys.exit(1)

# 2. PHYSICAL CONSTANTS
SIGMA_STRUCT = 11.5  
SIGMA_INT = 3.0      
THRESHOLD_SI = 287.35 

def load_sparc_mrt(path, columns, skip_lines):
    if not os.path.exists(path): raise FileNotFoundError(f"Missing: {path}")
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln for ln in f.readlines()[skip_lines:] if ln.strip() and not ln.strip().startswith(('#', '|')) and not set(ln.strip()) <= set("-=._ ")]
    df = pd.read_csv(io.StringIO("".join(lines)), sep=r"\s+", names=columns, engine="python", na_values=["----", "...", "NaN", "nan", ""])
    id_col = 'Galaxy' if 'Galaxy' in df.columns else 'ID'
    df['Galaxy'] = df[id_col].astype(object).apply(lambda x: str(x).strip() if pd.notnull(x) else np.nan)
    df = df[df['Galaxy'] != 'nan'].dropna(subset=['Galaxy']).copy()
    for c in df.columns:
        if c not in ["Galaxy", "ID", "Ref", "f_D"]: df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

# 3. AUDIT EXECUTION
if __name__ == "__main__":
    DATA_DIR = os.path.join(BASE_DIR, "data")
    meta = load_sparc_mrt(os.path.join(DATA_DIR, "sparc_metadata.csv"), ['Galaxy','T','D','e_D','f_D','Inc','e_Inc','L36','e_L36','Reff','SBeff','Rdisk','SBdisk','MHI','RHI','Vflat','e_Vflat','Q','Ref'], 34)
    mm = load_sparc_mrt(os.path.join(DATA_DIR, "sparc_mass_models.csv"), ['ID','D_mm','R','Vobs','errV','Vgas','Vdisk','Vbul','SBdisk','SBbul'], 25)

    df_merged = mm.merge(meta, on="Galaxy", how="inner")
    df = df_merged[(df_merged["Q"] == 1) & (df_merged["Inc"] >= 40.0) & (df_merged["Inc"] <= 70.0)].copy()
    df = df.groupby("Galaxy").filter(lambda g: len(g) >= 8).copy()

    valid_p, valid_t, valid_w = [], [], []
    phase_counts = {"Locked": 0, "Fluid": 0}
    sigma_vals = []

    print(f"📊 Auditing {df['Galaxy'].nunique()} Galaxies (Strict Blind Protocol)...")

    for gal, group in df.groupby("Galaxy"):
        # Protocol: Sort radii descending for Outside-In Resolution
        group = group.sort_values("R", ascending=False)
        R, Vobs, errV = group["R"].to_numpy(), group["Vobs"].to_numpy(), group["errV"].to_numpy()
        Vgas, Vdisk, Vbul = group["Vgas"].to_numpy(), group["Vdisk"].to_numpy(), np.nan_to_num(group["Vbul"].to_numpy())
        L36, MHI, Reff, RHI = group["L36"].iloc[0], group["MHI"].iloc[0], max(group["Reff"].iloc[0], 0.1), max(group["RHI"].iloc[0], 0.1)

        n_split = max(int(len(R) * 0.25), 2)
        w_pts = 1.0 / (errV**2 + SIGMA_INT**2 + SIGMA_STRUCT**2)

        # 1. ANCHORING (Outer 25%) - REMOVED GLOBAL LEAKAGE
        best_ud, best_t_boundary, best_chi2 = 0.5, 0.025, np.inf
        for ud_trial in [0.5, 0.7]:
            # Sigma calculation (Distance Invariant)
            trial_sigma = (ud_trial * L36 * 1e9 / (np.pi * (Reff * 1000)**2)) + (1.33 * MHI * 1e9 / (np.pi * (RHI * 1000)**2))
            logI_comp = np.clip(-0.05 + 0.20 * np.log10(max(trial_sigma, 1e-3)), 0.25, 0.60)
            
            # Newtonian velocity baseline for anchor set
            v_bar_sq_anchor = np.maximum(Vgas[:n_split]**2 + ud_trial*Vdisk[:n_split]**2 + 1.4*ud_trial*Vbul[:n_split]**2, 1e-12)
            # Use weighted mean of ONLY anchor points for initialization
            v_obs_anchor_mean = np.average(Vobs[:n_split], weights=w_pts[:n_split])
            
            local_t, min_err = 0.025, np.inf
            # Calibration loop using ONLY anchor-set statistics
            for t_trial in np.linspace(0.001, 0.45, 150):
                # Predict mean velocity for anchor zone
                vp = engine.predict_from_invariants(float(np.mean(np.sqrt(v_bar_sq_anchor))), float(t_trial), logI_comp)
                if abs(vp - v_obs_anchor_mean) < min_err: 
                    min_err, local_t = abs(vp - v_obs_anchor_mean), t_trial
            
            # Verify anchor fit
            vp_a = [engine.predict_from_invariants(float(np.sqrt(v_bar_sq_anchor[j])), float(local_t), logI_comp) for j in range(n_split)]
            chi = np.sum(w_pts[:n_split] * (Vobs[:n_split] - np.array(vp_a))**2)
            if chi < best_chi2: 
                best_chi2, best_ud, best_t_boundary = chi, ud_trial, local_t

        # 2. STATE REGISTRATION (SI Thresholding)
        final_sigma = (best_ud * L36 * 1e9 / (np.pi * (Reff * 1000)**2)) + (1.33 * MHI * 1e9 / (np.pi * (RHI * 1000)**2))
        sigma_vals.append(final_sigma)
        is_locked = final_sigma > THRESHOLD_SI
        k_state = 0.005 if is_locked else 0.012
        phase_counts["Locked" if is_locked else "Fluid"] += 1
        logI_comp = np.clip(-0.05 + 0.20 * np.log10(max(final_sigma, 1e-3)), 0.25, 0.60)

        # 3. BLIND PREDICTION (Inner 75%)
        v_bar = np.sqrt(np.maximum(Vgas**2 + best_ud*Vdisk**2 + 1.4*best_ud*Vbul**2, 1e-12))
        R0 = max(R[0], 0.05) 
        accel = v_bar**2 / np.maximum(R, 0.01) 
        
        for j in range(len(R)):
            # Hierarchy reduction logic applied blindly
            alpha = 1.0 - (R[j] / R0)
            red = k_state * alpha * (np.log10(max(accel[j] / accel[0], 1.0))**2)
            vp = engine.predict_from_invariants(float(v_bar[j]), float(max(best_t_boundary - red, 0.001)), logI_comp)
            
            # Record ONLY validation (blind) points for final metrics
            if j >= n_split:
                valid_p.append(vp); valid_t.append(Vobs[j]); valid_w.append(w_pts[j])

    # 4. SCIENTIFIC REPORTING
    y_p, y_t, w = np.array(valid_p), np.array(valid_t), np.array(valid_w)
    
    def weighted_corr(x, y, weights):
        avg_x = np.average(x, weights=weights)
        avg_y = np.average(y, weights=weights)
        diff_x = x - avg_x
        diff_y = y - avg_y
        cov = np.sum(weights * diff_x * diff_y)
        var_x = np.sum(weights * diff_x**2)
        var_y = np.sum(weights * diff_y**2)
        return cov / np.sqrt(var_x * var_y)

    r_pearson = weighted_corr(y_p, y_t, w)
    r2_val = 1 - (np.sum(w * (y_t - y_p)**2) / np.sum(w * (y_t - np.average(y_t, weights=w))**2))
    chi2_n = np.sum(w * (y_t - y_p)**2) / len(y_t)

    print(f"\n🏆 FINAL AUDIT V13.7.2 (STRICT BLIND)")
    print(f"{'='*60}")
    print(f"Pearson r (Blind):   {r_pearson:.4f}")
    print(f"VALIDATION-ONLY R²:  {r2_val:.4f}")
    print(f"χ²/N (Blind Zone):   {chi2_n:.4f}")
    print(f"{'='*60}")