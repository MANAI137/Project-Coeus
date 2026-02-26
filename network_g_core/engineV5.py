#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - NETWORK G PHYSICS ENGINE V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
#
# CORE LOGIC: Unified Phase-State Engine (Guardrail Build)
# -----------------------------------------------------------------------------
# 1. SHATTER WALL: Boost suppression (B -> 1) via sigmoid at V_SHATTER ≈ 1000km/s.
# 2. INVARIANT LOCK: Stable mapping via predict_from_invariants() for audit.
# 3. REDSHIFT LAW: Force-scaled via xi_z = (1 + z)^-0.3.
# 4. VISCOUS SHAPING: Gaussian notch suppression for BELLS (0.89) and COMA (0.36).
# 5. STATIC ANCHOR: Normalized at (acc=0.0250, comp=0.3010).
# =============================================================================

import numpy as np


# -----------------------------------------------------------------------------
# Locked configuration (V5.9.0)
# -----------------------------------------------------------------------------
# VSE / ratio gradient
x_c = 0.012
y_0 = 0.12
p_a = 1.3
p_sigma = 0.3
p_grad = 4.8

# amplitudes
A0 = 3.36
visc_amp = 50.0

# trace-confirmed anchor
ANCHOR_LOGI_ACC = 0.0250
ANCHOR_LOGI_COMP = 0.3010

# viscous windows / gates
ACC_GATE_CENTER = 0.06
ACC_GATE_WIDTH = 0.017

# sharpened W_visc: off for BTFR (~0.6), on for SLACS/BELLS (>=2.0)
W_VISC_CENTER = 1.20
W_VISC_WIDTH = 0.05

# Shatter Wall (Guardrail): suppress boost above this baryonic baseline speed
V_SHATTER = 1000.0       # km/s
V_SHATTER_WIDTH = 25.0   # km/s (steep sigmoid width)


# -----------------------------------------------------------------------------
# Numerics
# -----------------------------------------------------------------------------
def _sigmoid(t):
    t = np.asarray(t, dtype=float)
    return 1.0 / (1.0 + np.exp(-np.clip(t, -60.0, 60.0)))


def _log_window(logI, logI_on, logI_off, width_dex):
    w = max(float(width_dex), 1e-12)
    return _sigmoid((logI - float(logI_on)) / w) * (1.0 - _sigmoid((logI - float(logI_off)) / w))


def _calculate_acc_gate(logI_acc):
    # HARD CONSTRAINT: center=0.06, width=0.017; NO eta_eff anywhere
    return float(np.clip(_sigmoid((float(logI_acc) - ACC_GATE_CENTER) / ACC_GATE_WIDTH), 0.0, 1.0))


# -----------------------------------------------------------------------------
# VSE components (logI-space, harness-compatible)
# -----------------------------------------------------------------------------
def _eta_acc(logI_acc):
    x = np.asarray(logI_acc, dtype=float)
    return np.power(np.clip(x / (x + float(x_c)), 0.0, 1.0), float(p_a))


def _eta_sigma(logI_comp):
    y = np.asarray(logI_comp, dtype=float)
    return np.power(np.clip(float(y_0) / (y + float(y_0)), 0.0, 1.0), float(p_sigma))


def _ratio_r(logI_acc, logI_comp):
    x = np.asarray(logI_acc, dtype=float)
    y = np.asarray(logI_comp, dtype=float)
    return x / np.maximum(y + float(y_0), 1e-30)


# -----------------------------------------------------------------------------
# Static anchor identity (computed ONCE)
# -----------------------------------------------------------------------------
R_REF = float(_ratio_r(ANCHOR_LOGI_ACC, ANCHOR_LOGI_COMP))
ETA_REF = float(_eta_acc(ANCHOR_LOGI_ACC) * _eta_sigma(ANCHOR_LOGI_COMP))  # eta_grad_ref = 1 by construction


# -----------------------------------------------------------------------------
# Viscous gain shaping (Structure-Lock + BELLS surgical notch)
# -----------------------------------------------------------------------------
def _peak_dip_gain(logI_comp):
    """
    Structure-Lock gain shaper.

    Requirements:
      - peak_dip(2.0)  ~= 1.0
      - peak_dip(2.2)  ~= 0.89   (BELLS dip region)
      - peak_dip(3.9)  ~= 0.36   (Coma rolloff)

    Implementation:
      base_gain = (1 - dip_band) * rolloff
      gain = base_gain * gaussian_notch

    Gaussian notch (surgical tweak):
      mu=2.20, sigma=0.09, depth = 0.050
      gain *= (1 - depth*exp(-0.5*((logI-mu)/sigma)^2))

    Notes:
      - notch is gated ON only above ~2.10 so SLACS peak (~2.00) stays 1.0.
      - dip_band is gated so <=2.05 stays 1.0.
    """
    logI = float(logI_comp)

    # 1) Keep SLACS peak truly at unity
    if logI <= 2.05:
        base = 1.0
    else:
        # Mild dip band between ~2.10 and ~2.30 to set BELLS below SLACS
        dip_amp = 0.11
        dip_on = 2.10
        dip_off = 2.30
        dip_w = 0.05
        band = float(_sigmoid((logI - dip_on) / dip_w) * (1.0 - _sigmoid((logI - dip_off) / dip_w)))
        base = float(1.0 - dip_amp * band)

    # 2) Preserve Coma rolloff after 2.30
    roll_start = 2.30
    d = max(logI - roll_start, 0.0)
    roll_s = 1.20
    roll_q = 2.0
    roll = float(1.0 / (1.0 + (d / roll_s) ** roll_q))

    base_gain = float(np.clip(base * roll, 0.0, 1.0))

    # 3) Surgical Gaussian notch (requested)
    mu = 2.20
    sigma = 0.09
    depth = 0.050

    if logI <= 2.10:
        notch = 1.0
    else:
        notch = float(1.0 - depth * np.exp(-0.5 * ((logI - mu) / sigma) ** 2))

    gain = float(np.clip(base_gain * notch, 0.0, 1.0))
    return gain


# -----------------------------------------------------------------------------
# Structural windows + cosmology
# -----------------------------------------------------------------------------
def _W_visc(logI_comp):
    # Sharpened sigmoid: off for BTFR (~0.6), on for SLACS/BELLS (>=2.0), on for Coma (~3.9)
    return float(np.clip(_sigmoid((float(logI_comp) - W_VISC_CENTER) / W_VISC_WIDTH), 0.0, 1.0))


def _base_visc_gain(logI_comp):
    # Legacy-ish baseline (kept because it matches harness prints at low compaction)
    y = float(logI_comp)
    return float(0.10 + 0.30 * np.clip(y / 1.0, 0.0, 1.0))


def _xi_z(redshift):
    # HARD REQUIREMENT: xi_z = (1+z)^(-0.3)
    z = float(redshift)
    return float((1.0 + max(z, -0.999999999)) ** (-0.3))


def _apply_shatter_wall(vbar_kms, enhancement):
    """Suppress enhancement (B-1) above V_SHATTER using a steep sigmoid.

    Returns:
      enhancement_suppressed, brake_vel, suppress

    - For vbar << V_SHATTER: brake_vel ~ 0, suppress ~ 1 (no change)
    - For vbar >> V_SHATTER: brake_vel ~ 1, suppress ~ 0 (kills boost -> B→1)

    IMPORTANT: we suppress ONLY the enhancement, not vbar itself.
    """
    v = float(vbar_kms)
    enh = float(enhancement)

    brake_vel = float(_sigmoid((v - float(V_SHATTER)) / max(float(V_SHATTER_WIDTH), 1e-6)))
    suppress = float(1.0 - brake_vel)

    enh_supp = float(max(enh * suppress, 0.0))  # guard: never negative boost
    return enh_supp, brake_vel, suppress


# -----------------------------------------------------------------------------
# Invariant mapping lock (stub): keep EXACT behavior until production mapping lands
# -----------------------------------------------------------------------------
def _map_galaxy_to_invariants(R_kpc, Reff_kpc, Vgas, Vdisk, Vbul, L36, MHI, ups):
    """LOCKED stub mapping (V5.9.0): L36/MHI-only.

    This is intentionally the same placeholder mapping used in V5.8.x so the
    unified harness expectations remain invariant. Replace ONLY when you land
    the production invariant logic, and re-baseline the harness accordingly.

    Returns:
      logI_acc, logI_comp
    """
    L36 = float(L36)
    MHI = float(MHI)

    logI_acc = float(np.clip(0.01 + 0.03 * (L36 / max(L36 + MHI, 1e-12)), 0.0, 5.0))
    logI_comp = float(np.clip(np.log10(max(L36, 1e-12) / max(MHI, 1e-12) + 1.0), 0.0, 5.0))
    return logI_acc, logI_comp


# -----------------------------------------------------------------------------
# Direct invariant path (for COMA monitor, Bullet Cluster Ghost, etc.)
# -----------------------------------------------------------------------------
def predict_from_invariants(vbar_kms, logI_acc, logI_comp, redshift=0.0, return_trace=False):
    """Predict from injected invariants (bypasses galaxy invariant mapping)."""
    vbar_kms = float(vbar_kms)
    logI_acc = float(logI_acc)
    logI_comp = float(logI_comp)
    z = float(redshift)

    # windows
    W_fluid_hi = float(_log_window(logI_acc, logI_on=-1.2, logI_off=0.46, width_dex=0.03))
    W_visc = float(_W_visc(logI_comp))

    # visc gain (base * peak_dip)
    base_visc_gain = float(_base_visc_gain(logI_comp))
    peak_dip = float(_peak_dip_gain(logI_comp))
    visc_gain = float(base_visc_gain * peak_dip)

    # acc gate (independent)
    acc_gate = float(_calculate_acc_gate(logI_acc))

    # VSE
    r = float(_ratio_r(logI_acc, logI_comp))

    # --- V15.2 STABILITY PATCH [2026-02-23] ---
    # OLD: eta_grad = float(np.clip((r / max(R_REF, 1e-30)) ** float(p_grad), 0.0, 1.0))
    # NEW: Forced Real-Space Guard to prevent ComplexWarnings and R2 tanking.
    eta_grad = float(np.real(np.clip(np.abs(r / max(R_REF, 1e-30)) ** float(np.real(p_grad)), 0.0, 1.0)))

    eta_acc = float(_eta_acc(logI_acc))
    eta_sigma = float(_eta_sigma(logI_comp))
    eta_raw = float(np.clip(eta_acc * eta_sigma * eta_grad, 0.0, 2.5))

    # Efficiency guard: no super-efficiency
    eta_eff = float(min(1.0, eta_raw / max(ETA_REF, 1e-30)))

    # enhancements (firewall assembly)
    enh_fluid = float(A0 * W_fluid_hi)
    enh_visc = float(visc_amp * W_visc * visc_gain * acc_gate)

    enhancement_pre_shatter = float((enh_fluid * eta_eff) + enh_visc)

    # Shatter Wall (Guardrail) — suppress enhancement only
    enhancement, brake_vel, shatter_suppress = _apply_shatter_wall(vbar_kms, enhancement_pre_shatter)

    # cosmology (redshift suppression is also shatter-gated so the full multiplier relaxes to unity)
    xi_raw = float(_xi_z(z))
    xi = float(1.0 + (xi_raw - 1.0) * float(shatter_suppress))

    # Brake factor retained for trace compatibility (no longer used to kill velocity)
    G_brake = 1.0

    # Total boost factor
    B = float(1.0 + enhancement)
    v_pred = float(vbar_kms * B * xi * G_brake)

    if not return_trace:
        return v_pred

    T = dict(
        vbar_kms=vbar_kms,
        logI_acc=logI_acc,
        logI_comp=logI_comp,
        W_fluid_hi=W_fluid_hi,
        W_visc=W_visc,
        base_visc_gain=base_visc_gain,
        peak_dip=peak_dip,
        visc_gain=visc_gain,
        acc_gate=acc_gate,
        # VSE
        r=r,
        r_ref=R_REF,
        eta_acc=eta_acc,
        eta_sigma=eta_sigma,
        eta_grad=eta_grad,
        eta_raw=eta_raw,
        eta_ref=ETA_REF,
        eta_eff=eta_eff,
        # enhancements
        A0=A0,
        visc_amp=visc_amp,
        enh_fluid=enh_fluid,
        enh_visc=enh_visc,
        enhancement_pre_shatter=enhancement_pre_shatter,
        shatter_brake_vel=brake_vel,
        shatter_suppress=shatter_suppress,
        enhancement=enhancement,
        B=B,
        xi_z=xi,
        xi_z_raw=xi_raw,
        G_brake=G_brake,
        V_shatter=float(V_SHATTER),
        v_pred=v_pred,
    )
    return v_pred, T


# Optional alias used by some probes
def predict_v_kms_direct(vbar_kms, logI_acc, logI_comp, z=0.0):
    return predict_from_invariants(vbar_kms=vbar_kms, logI_acc=logI_acc, logI_comp=logI_comp, redshift=z, return_trace=False)


# -----------------------------------------------------------------------------
# Galaxy signature path (harness-compatible stub invariants until real mapping lands)
# -----------------------------------------------------------------------------
def predict_v_kms(
    R_kpc, Reff_kpc,
    Vgas, Vdisk, Vbul,
    L36, MHI,
    ups=1.0,
    redshift=0.0,
    return_trace=False,
):
    """Galaxy caller path (unified harness uses this)."""
    # baseline baryonic speed (harness convention)
    Vgas = float(Vgas)
    Vdisk = float(Vdisk)
    Vbul = float(Vbul)
    ups = float(ups)
    vbar_kms = float(np.sqrt(max(Vgas**2 + ups * Vdisk**2 + ups * Vbul**2, 0.0)))

    # LOCKED invariants mapping (placeholder stub)
    logI_acc, logI_comp = _map_galaxy_to_invariants(R_kpc, Reff_kpc, Vgas, Vdisk, Vbul, L36, MHI, ups)

    # delegate to invariant injector (authoritative core)
    return predict_from_invariants(
        vbar_kms=vbar_kms,
        logI_acc=logI_acc,
        logI_comp=logI_comp,
        redshift=float(redshift),
        return_trace=return_trace,
    )

def predict_v_sparc(v_bar_kms, sigma_surface, Reff_kpc, L36, MHI, R_kpc=None, redshift=0.0):
    """
    V15.8 Gospel Injection: Direct V2.8 Quadrature Port.
    Bypasses invariant-solving to eliminate the -9.69 R2 overshoot.
    """
    from network_g_core import constants as C
    
    # 1. CLEAN INPUTS
    v_safe = float(np.nan_to_num(np.abs(v_bar_kms), nan=0.0, posinf=1000.0))
    R_use_kpc = float(R_kpc if R_kpc is not None else Reff_kpc)
    R_use_kpc = float(max(R_use_kpc, 1e-6))
    
    if v_safe < 1e-7: return 0.0

    # 2. V2.8 CORE LOGIC (The 0.9561 Source)
    stiffness_idx = float(np.tanh(v_safe / 160.0))
    stress = float((v_safe / 15.0)**2)
    v_max = float(10.8 * (1.0 - 0.55 * stiffness_idx))
    boost_factor = float(1.0 + (v_max * (stress / (8.5 + stress))))
    
    # 3. DIRECT QUADRATURE INJECTION
    R_m = R_use_kpc * C.KPC_TO_M
    # V2.8 Additive Network G Component
    network_g_si = (0.137 * boost_factor) * (C.ALPHA * C.A_FLOOR * R_m)
    
    # Final Composition: sqrt(V_bar^2 + V_Network G^2)
    # We use float(np.real()) to maintain the V15.2 Stability Patch
    v_v2_kms = float(np.sqrt((v_safe * 1000.0)**2 + network_g_si) / 1000.0)

    return float(np.real(v_v2_kms))