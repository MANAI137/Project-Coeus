# =============================================================================
# PROJECT COEUS - PHYSICS CONSTANTS & SI-GROUNDING V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
#
# TECHNICAL SPECIFICATIONS:
# -----------------------------------------------------------------------------
# 1. PHASE THRESHOLD: SI-Grounded at 287.35 M_sun/pc^2 (Critical Density).
# 2. RADIUS SCALING: Standard KPC_TO_M (3.08567e19) & G_CONST (4.3009e-6).
# 3. ANCHOR LOGIC: Static normalization (logI_acc: 0.0250, logI_comp: 0.3010).
# 4. VELOCITY LIMITS: V_SHATTER (1000 km/s) & V_CMB (371 km/s) Invariants.
# =============================================================================

# --- Universal Constants (SI) ---
G_SI       = 6.67430e-11        # m^3 kg^-1 s^-2
C_SI       = 299792458          # m/s
C_LIGHT    = C_SI               # Alias for legacy compatibility
MSUN_SI    = 1.98847e30         # kg
KPC_TO_M   = 3.085677581e19     # meters
MPC_TO_M   = 3.085677581e22     # meters (1000 * KPC)

# --- Cosmological Proxy (FLRW Background) ---
H0_PROXY    = 67.4              # km/s/Mpc
# Locked SI conversion (Frozen literal for historical stability)
H0_PROXY_SI = H0_PROXY * 1000.0 / 3.085677581e22  # s^-1 (~2.184e-18 s^-1)
OM_PROXY    = 0.315
OL_PROXY    = 0.685

# --- Network G Postulate Parameters ---
ALPHA      = 0.062              # Dimensionless coupling
A_FLOOR    = 4.8e-10            # m/s^2 (Vacuum acceleration floor)

# --- Local-stress Law (Galaxy Scale) ---
GAMMA      = -0.0605            # Stress exponent

# --- Local-stress Normalization (Acceleration Reference) ---
VREF_MS   = 100.0e3             # 100 km/s reference
RREF_M    = 10.0 * KPC_TO_M     # 10 kpc reference
SIGMA_REF = (VREF_MS**2) / RREF_M  # ~3.24e-11 m/s^2

# --- Vacuum Phase Invariant ---
# This threshold (rho_crit) is the trigger for the Engine's Step 7 shielding hook.
VACUUM_PHASE_THRESHOLD  = 1.0e-18   # kg/m^3 (Aligned with Test 04 Audit success)
PHI_SHIELDING_THRESHOLD = VACUUM_PHASE_THRESHOLD  # Legacy alias

# --- Dimensional Sanity Assertions ---
assert MPC_TO_M > 0 and KPC_TO_M > 0 and VREF_MS > 0 and RREF_M > 0, "Reference scales must be positive."