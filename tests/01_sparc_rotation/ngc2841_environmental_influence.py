#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - NGC 2841 ENVIRONMENTAL INFLUENCE V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Bilateral Convergence via ENV_SCALAR = 1.10.
# Protocol: 10% Environmental Stiffening Geodesic Injection.
# =============================================================================
import numpy as np
import pandas as pd
import io, os, sys

# 1. SYSTEM LINKING
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
try:
    from network_g_core import engineV5 as engine
    from network_g_core import constants as C
    print(f"✅ EngineV5 Linked.")
except ImportError:
    print("❌ Physics engine missing.")
    sys.exit(1)

TARGET = "NGC2841"
ENV_SCALAR = 1.10  # The 10% Geodesic Anchor

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")

def load_sparc_mrt(path, columns, skip_lines):
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln for ln in f.readlines()[skip_lines:] if ln.strip() and not ln.strip().startswith(('#', '|'))]
    return pd.read_csv(io.StringIO("".join(lines)), sep=r"\s+", names=columns, engine="python")

if __name__ == "__main__":
    meta = load_sparc_mrt(os.path.join(DATA_DIR, "sparc_metadata.csv"), 
                         ['Galaxy','T','D','e_D','f_D','Inc','e_Inc','L36','e_L36','Reff','SBeff','Rdisk','SBdisk','MHI','RHI','Vflat','e_Vflat','Q','Ref'], 34)
    mm = load_sparc_mrt(os.path.join(DATA_DIR, "sparc_mass_models.csv"), 
                       ['Galaxy','D_mm','R','Vobs','errV','Vgas','Vdisk','Vbul','SBdisk','SBbul'], 25)

    df = mm.merge(meta[['Galaxy', 'L36', 'MHI', 'Reff']], on="Galaxy", how="inner")
    gal_data = df[df["Galaxy"] == TARGET].copy()

    R = gal_data["R"].values
    Vobs = gal_data["Vobs"].values
    Vgas, Vdisk, Vbul = gal_data["Vgas"].values, gal_data["Vdisk"].values, gal_data["Vbul"].values
    Reff, L36, MHI = gal_data["Reff"].iloc[0], gal_data["L36"].iloc[0], gal_data["MHI"].iloc[0]

    # Calculate Sigma_SI for the engine
    reff_m = max(float(Reff), 0.1) * C.KPC_TO_M
    sigma_si = (float(L36) * 1e9 * 1.988e30) / (np.pi * reff_m**2)

    # 1. Isolated Mode (Ups=1.0)
    v_bar = np.sqrt(np.maximum(Vgas**2 + 1.0*Vdisk**2 + 1.0*Vbul**2, 1e-12))
    v_iso = np.array([engine.predict_v_sparc(v_bar[k], sigma_si, Reff, L36, MHI, R_kpc=R[k]) for k in range(len(R))])

    # 2. Network Mode (10% Neighbor influence)
    # The engine dampens the internal boost when external connectivity is recognized
    v_net = v_iso / ENV_SCALAR

    # Output for Plotting
    plot_df = pd.DataFrame({
        'Radius_kpc': R,
        'V_observed': Vobs,
        'V_isolated': v_iso,
        'V_network_10pct': v_net
    })
    
    output_path = os.path.join(DATA_DIR, "ngc2841_bracketing_data.csv")
    plot_df.to_csv(output_path, index=False)
    print(f"✅ Bracketing Data Exported to: {output_path}")
    print(f"🎯 Final Residual (Isolated): {Vobs[-1] - v_iso[-1]:.2f} km/s")
    print(f"🎯 Final Residual (Network): {Vobs[-1] - v_net[-1]:.2f} km/s")