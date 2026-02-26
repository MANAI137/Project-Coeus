# =============================================================================
# PROJECT COEUS - UNIFIED LAW VALIDATION HARNESS V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: 4-Point Regime Trace-Based Validation (SPARC/SLACS/BELLS/BULLET).
# Protocol: Multi-Regime Sync & COMA Cluster Invariant Injection.
# =============================================================================

import sys
import numpy as np

import network_g_core.engineV5 as engine
from network_g_core import constants

ALLOWED_KEYS = {"R_kpc", "Reff_kpc", "Vgas", "Vdisk", "Vbul", "L36", "MHI", "ups", "redshift"}


# -----------------------------------------------------------------------------
# Trace capture: prefer engine-native trace, else sys profiler fallback
# -----------------------------------------------------------------------------
def capture_predict_trace(**kwargs):
    """
    Returns:
      v_kms (float), trace (dict)

    Works with:
      - V5.+ engines that implement return_trace=True
    """
    engine_params = {k: kwargs[k] for k in ALLOWED_KEYS if k in kwargs}

    # Preferred: engine-native trace
    try:
        out = engine.predict_v_kms(**engine_params, return_trace=True)
        if isinstance(out, (tuple, list)) and len(out) == 2 and isinstance(out[1], dict):
            return float(out[0]), dict(out[1])
    except TypeError:
        pass
    except Exception:
        pass

    # Fallback: sys profiler capture of locals
    target_code = engine.predict_v_kms.__code__
    captured = {"locals": None}

    def profiler(frame, event, arg):
        if event == "return" and frame.f_code is target_code:
            captured["locals"] = dict(frame.f_locals)
        return profiler

    old = sys.getprofile()
    sys.setprofile(profiler)
    try:
        v = engine.predict_v_kms(**engine_params)
    finally:
        sys.setprofile(old)

    return float(v), (captured["locals"] or {})


def compute_vbar_kms(Vgas, Vdisk, Vbul, ups):
    # NOTE: vbar is the kinematic scalar baseline in km/s (not squared, not SI)
    return float(np.sqrt(max(float(Vgas) ** 2 + float(ups) * float(Vdisk) ** 2 + float(ups) * float(Vbul) ** 2, 0.0)))


# Precision map for diagnostics
_PREC = {
    "eta_eff": 6,
    "visc_gain": 6,
    "shatter_brake_vel": 6,
    "shatter_suppress": 6,
}


def _fmt_val(d, key, w=9, prec=2):
    if key not in d:
        return f"{'NA':>{w}}"
    try:
        p = _PREC.get(key, prec)
        return f"{float(d[key]):{w}.{p}f}"
    except Exception:
        s = str(d[key])
        return f"{s[:w]:>{w}}"


def _pick(d, *keys):
    """Return the first present key in d, else None."""
    for k in keys:
        if k in d:
            return k
    return None


def run_case(name, params, target_boost=None, target_v=None, boost_ref="vbar"):
    """
    boost_ref:
      - "vbar"  compares target_boost to boost(v/vbar)
      - "Vdisk" compares target_boost to boost(v/Vdisk_input)
    """
    engine_params = {k: params[k] for k in ALLOWED_KEYS if k in params}

    v_engine, T = capture_predict_trace(**engine_params)

    Vgas = float(engine_params["Vgas"])
    Vdisk = float(engine_params["Vdisk"])
    Vbul = float(engine_params["Vbul"])
    ups = float(engine_params["ups"])

    vbar_kms = compute_vbar_kms(Vgas, Vdisk, Vbul, ups)

    boost_vs_vbar = v_engine / max(vbar_kms, 1e-12)
    boost_vs_Vdisk = v_engine / max(Vdisk, 1e-12)

    boost_for_req = boost_vs_vbar if boost_ref == "vbar" else boost_vs_Vdisk

    print(name)
    if target_boost is not None:
        print(f"  target boost ({boost_ref}): {target_boost:6.2f} | engine boost: {boost_for_req:6.2f}")
    if target_v is not None:
        print(f"  target v:            {target_v:8.2f} | engine v: {v_engine:8.2f}")

    print(
        f"  [BASE] Vdisk_in: {Vdisk:8.2f} | vbar_kms: {vbar_kms:8.2f} | "
        f"boost(v/vbar): {boost_vs_vbar:6.2f} | boost(v/Vdisk): {boost_vs_Vdisk:6.2f}"
    )

    # Diagnostics: common fields
    k_Wf = _pick(T, "W_fluid", "W_fluid_hi")
    Wf_label = "W_fluid" if k_Wf == "W_fluid" else "W_fluid_hi"
    Wf_val = _fmt_val(T, k_Wf, w=7) if k_Wf else f"{'NA':>7}"

    print(
        "  [DIAG]"
        f" logI_acc:{_fmt_val(T,'logI_acc', w=9, prec=4)}"
        f" | logI_comp:{_fmt_val(T,'logI_comp', w=9, prec=4)}"
        f" | {Wf_label}:{Wf_val}"
        f" | W_visc:{_fmt_val(T,'W_visc', w=9)}"
        f" | visc_gain:{_fmt_val(T,'visc_gain', w=11)}"
    )

    enh_key = _pick(T, "enhancement", "enh")
    enh_val = _fmt_val(T, enh_key, w=9) if enh_key else f"{'NA':>9}"

    print(
        "  [DIAG]"
        f" enh:{enh_val}"
        f" | B:{_fmt_val(T,'B', w=9)}"
        f" | xi_z:{_fmt_val(T,'xi_z', w=9)}"
        f" | G_brake:{_fmt_val(T,'G_brake', w=9)}"
    )

    # VSE diagnostics
    if any(k in T for k in ("eta_acc", "eta_sigma", "eta_grad", "eta_eff", "enh_fluid", "enh_visc")):
        print(
            "  [VSE ]"
            f" eta_acc:{_fmt_val(T,'eta_acc', w=9)}"
            f" | eta_sigma:{_fmt_val(T,'eta_sigma', w=9)}"
            f" | eta_grad:{_fmt_val(T,'eta_grad', w=9)}"
            f" | eta_eff:{_fmt_val(T,'eta_eff', w=11)}"
        )
        print(
            "  [VSE ]"
            f" enh_fluid:{_fmt_val(T,'enh_fluid', w=9)}"
            f" | enh_visc:{_fmt_val(T,'enh_visc', w=9)}"
            f" | acc_gate:{_fmt_val(T,'acc_gate', w=9)}"
            f" | visc_amp:{_fmt_val(T,'visc_amp', w=9)}"
        )

    # Shatter diagnostics (if present)
    if any(k in T for k in ("shatter_brake_vel", "shatter_suppress", "V_shatter")):
        print(
            "  [SHAT]"
            f" V_shatter:{_fmt_val(T,'V_shatter', w=9)}"
            f" | brake_vel:{_fmt_val(T,'shatter_brake_vel', w=11)}"
            f" | suppress:{_fmt_val(T,'shatter_suppress', w=11)}"
        )

    print("-" * 94)


# -----------------------------------------------------------------------------
# 4-point regime harness
# -----------------------------------------------------------------------------
def run_4point_harness(mode="CANONICAL"):
    """
    mode:
      - "CANONICAL": Vdisk_in fixed at 100 for SPARC/SLACS/BELLS (legacy comparable).
      - "VBAR_NORMALIZED": scales Vdisk so vbar=100 for SPARC/SLACS/BELLS (changes invariants).
    """
    print("🧪 4-POINT REGIME HARNESS (TRACE-BASED, DUAL-BOOST)")
    print(f"Mode: {mode}")
    print("═" * 94)

    z = 0.2
    V_in = 100.0  # legacy harness "V" input

    def Vdisk_for(ups):
        if mode == "VBAR_NORMALIZED":
            return V_in / np.sqrt(ups)
        return V_in

    cases = [
        ("SPARC (Fluid)", dict(
            R_kpc=8.0, Reff_kpc=3.0,
            Vgas=0.0, Vbul=0.0,
            ups=1.0,
            Vdisk=Vdisk_for(1.0),
            L36=5.0, MHI=5.0,
            redshift=z
        ), 4.13),

        ("SLACS (Viscous)", dict(
            R_kpc=5.0, Reff_kpc=2.0,
            Vgas=0.0, Vbul=0.0,
            ups=1.2,
            Vdisk=Vdisk_for(1.2),
            L36=50.0, MHI=0.5,
            redshift=z
        ), 5.88),

        ("BELLS (Viscous)", dict(
            R_kpc=7.0, Reff_kpc=2.5,
            Vgas=0.0, Vbul=0.0,
            ups=1.2,
            Vdisk=Vdisk_for(1.2),
            L36=80.0, MHI=0.5,
            redshift=z
        ), 5.26),
    ]

    # Bullet: should trip Shatter Wall (boost ~ 1; enhancement ~ 0)
    V_bullet = 3000.0
    cases.append((
        "Bullet (Solid) — Shatter Wall",
        dict(
            R_kpc=300.0, Reff_kpc=250.0,
            Vgas=0.0, Vbul=0.0,
            ups=1.0,
            Vdisk=V_bullet,
            L36=50000.0, MHI=50000.0,
            redshift=z
        ),
        1.00
    ))

    boost_ref_for_req = "Vdisk"  # legacy tables used v_pred / 100
    for name, p, req in cases:
        run_case(name, p, target_boost=req, boost_ref=boost_ref_for_req)


# -----------------------------------------------------------------------------
# BTFR SWEEP (V12 dynamic-scaling, BASELINE-ALIGNED)
# -----------------------------------------------------------------------------
def run_btfr_sweep(size_scale=1.0):
    print("\n📊 BTFR SWEEP (V12 dynamic-scaling, BASELINE-ALIGNED)")
    print("═" * 94)
    print("Goal: Verify engine recovers BTFR slope with consistent baseline bookkeeping.")

    test_cases = [
        {"name": "Dwarf Spiral", "M": 1e8,  "V_target": 45},
        {"name": "Small Spiral", "M": 1e9,  "V_target": 80},
        {"name": "Mid Spiral",   "M": 1e10, "V_target": 145},
        {"name": "Giant Spiral", "M": 1e11, "V_target": 255},
    ]

    # Scaling Geometry (alpha=0.35 ensures surface density varies with mass)
    R0_kpc = 2.5
    alpha = 0.35  
    kappa = 2.5   
    f_enclosed = 0.60

    print(f"Assumptions: Reff = {R0_kpc}*(M/1e10)^{alpha}, Evaluation R = {kappa}*Reff")
    
    for t in test_cases:
        Mbar = float(t["M"])
        Vt = float(t["V_target"])

        # Reff grows slower than sqrt(M), so density (M/R^2) increases with Mass
        Reff_kpc = float(size_scale) * (R0_kpc * (Mbar / 1e10) ** alpha)
        R_kpc = kappa * Reff_kpc

        # Calculate Newtonian Vbar accurately
        Menc = f_enclosed * Mbar
        R_m = R_kpc * 3.086e19 # KPC_TO_M
        Vnewt_ms = np.sqrt(6.674e-11 * (Menc * 1.989e30) / np.maximum(R_m, 1e-30))
        Vbar_kms = Vnewt_ms / 1000.0

        # FIX: Ensure quadrature sum matches Vbar_kms (Removes baseline bias)
        # Partition: 30% Gas Velocity, ~95.4% Disk Velocity => sqrt(0.3^2 + 0.954^2) = 1.0
        Vgas = 0.30 * Vbar_kms
        Vdisk = 0.9539 * Vbar_kms

        # Run Engine Prediction
        v_pred, trace = capture_predict_trace(
            R_kpc=R_kpc, Reff_kpc=Reff_kpc,
            Vgas=Vgas, Vdisk=Vdisk, Vbul=0.0,
            L36=(0.6*Mbar/1e9), MHI=(0.4*Mbar/1e9),
            ups=1.0, redshift=0.02
        )

        boost = v_pred / max(Vbar_kms, 1e-10)
        err = (v_pred - Vt) / Vt * 100

        print(f"{t['name']:<13} | M: {Mbar:.1e} | vbar: {Vbar_kms:6.2f} | vpred: {v_pred:6.2f} | target: {Vt:3} | boost: {boost:.2f} | err: {err:+6.1f}%")
        
        if trace:
            # logI_comp will now vary from ~0.3 to ~1.4 across this sweep
            print(f"      [DIAG] logI_acc: {trace.get('logI_acc',0):.4f} | logI_comp: {trace.get('logI_comp',0):.4f} | Gate: {trace.get('W_fluid_hi',0):.1f}")


# -----------------------------------------------------------------------------
# COMA cluster monitor (invariant injection)
# -----------------------------------------------------------------------------
def run_coma_monitor():
    print("🔭 COMA CLUSTER MONITOR (INVARIANT-INJECTION)")
    print("═" * 94)

    # Physical baseline (from your harness)
    vbar_kms = 567.37
    z = 0.023

    # Injected invariants
    logI_acc = 0.038
    logI_comp = 3.92

    # Target: total virial velocity
    v_target = 1441.12
    boost_target = v_target / vbar_kms

    # Prefer engine's invariant injector if present
    try:
        out = engine.predict_from_invariants(
            vbar_kms=vbar_kms, logI_acc=logI_acc, logI_comp=logI_comp, redshift=z, return_trace=True
        )
        v_pred, T = float(out[0]), dict(out[1])
    except Exception:
        # Fallback to older engines if needed
        v_pred = float(engine.predict_v_kms_direct(vbar_kms, logI_acc, logI_comp, z=z))
        T = {}

    boost = v_pred / max(vbar_kms, 1e-12)

    print("COMA (Cluster)")
    print(f"  target boost (vbar): {boost_target:6.2f} | engine boost: {boost:6.2f}")
    print(f"  target v:            {v_target:8.2f} | engine v: {v_pred:8.2f}")
    print(f"  [BASE] vbar_kms: {vbar_kms:8.2f} | boost(v/vbar): {boost:6.2f}")
    print(
        "  [DIAG]"
        f" logI_acc:{logI_acc:9.4f}"
        f" | logI_comp:{logI_comp:9.4f}"
        f" | W_visc:{_fmt_val(T,'W_visc', w=9)}"
        f" | visc_gain:{_fmt_val(T,'visc_gain', w=11)}"
    )
    print(
        "  [DIAG]"
        f" enh:{_fmt_val(T,'enhancement', w=9)}"
        f" | B:{_fmt_val(T,'B', w=9)}"
        f" | xi_z:{_fmt_val(T,'xi_z', w=9)}"
        f" | G_brake:{_fmt_val(T,'G_brake', w=9)}"
    )
    if any(k in T for k in ("peak_dip", "base_visc_gain")):
        print(f"  [VISC] peak_dip:{_fmt_val(T,'peak_dip', w=9)} | base_visc_gain:{_fmt_val(T,'base_visc_gain', w=9)}")
    if any(k in T for k in ("eta_acc", "eta_sigma", "eta_grad", "eta_eff")):
        print(
            f"  [VSE ] eta_acc:{_fmt_val(T,'eta_acc', w=9)} | eta_sigma:{_fmt_val(T,'eta_sigma', w=9)}"
            f" | eta_grad:{_fmt_val(T,'eta_grad', w=9)} | eta_eff:{_fmt_val(T,'eta_eff', w=11)}"
        )

    print("-" * 94)


# -----------------------------------------------------------------------------
# TEST 06: Bullet Cluster Ghost (A vs B discrimination)
# -----------------------------------------------------------------------------
def run_bullet_ghost():
    print("🛰️  TEST 06: BULLET CLUSTER GHOST (A vs B)")
    print("═" * 94)

    # Keep vbar constant so discrimination is purely structural/invariant
    vbar = 100.0
    z = 0.2

    # Point A: Gas (High density / Low stress)
    A = dict(logI_acc=0.02, logI_comp=3.50)

    # Point B: Galaxies (Low density / High stress)
    B = dict(logI_acc=0.05, logI_comp=0.50)

    def pred_boost(pt):
        try:
            v = engine.predict_from_invariants(
                vbar_kms=vbar, logI_acc=pt["logI_acc"], logI_comp=pt["logI_comp"], redshift=z, return_trace=False
            )
        except Exception:
            v = engine.predict_v_kms_direct(vbar, pt["logI_acc"], pt["logI_comp"], z=z)
        return float(v) / vbar

    boost_A = pred_boost(A)
    boost_B = pred_boost(B)
    ratio = boost_B / max(boost_A, 1e-30)

    print(f"Point A (Gas):      logI_acc={A['logI_acc']:.4f}, logI_comp={A['logI_comp']:.4f} | boost={boost_A:.3f}")
    print(f"Point B (Galaxies): logI_acc={B['logI_acc']:.4f}, logI_comp={B['logI_comp']:.4f} | boost={boost_B:.3f}")
    print("─" * 78)
    print(f"Boost Ratio (B/A): {ratio:.3f}  -> {'✅ PASS' if ratio > 1.5 else '❌ FAIL'}")
    print("-" * 94)


def main():
    run_4point_harness(mode="CANONICAL")
    run_btfr_sweep(size_scale=1.0)
    run_coma_monitor()
    run_bullet_ghost()


if __name__ == "__main__":
    main()