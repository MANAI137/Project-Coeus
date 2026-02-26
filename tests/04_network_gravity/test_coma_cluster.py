#!/usr/bin/env python3
# =============================================================================
# PROJECT COEUS - COMA CLUSTER INVARIANT INJECTION V.1.0
# Author: Miguel Navarro
# ORCID: 0009-0009-5600-7985
# License: CC BY 4.0
# Logic: Viscous Transition in Isolated High-Mass Clusters.
# Protocol: Engine-Faithful Invariant Probe (test_coma_cluster.py).
# =============================================================================
import numpy as np
from network_g_core import engineV5  # engine module

def sigmoid(x):
    x = np.clip(x, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(-x))

def run_isolated_coma_probe():
    print("══════════════════════════════════════════════════════════════════════════════════════════")
    print("🧪 ISOLATED COMA PROBE (ENGINE-FAITHFUL INVARIANT INJECTION)")
    print("══════════════════════════════════════════════════════════════════════════════════════════")

    # --- Inputs (your chosen invariant point) ---
    v_bar = 567.0
    logI_acc  = 0.038
    logI_comp = 3.92
    z_coma = 0.023

    # --- Pull config safely (module attrs OR ENGINE_CFG dict) ---
    cfg = getattr(engineV5, "ENGINE_CFG", {})
    def pick(name, default=None):
        return float(getattr(engineV5, name, cfg.get(name, default)))

    A0       = pick("A0", 3.36)
    visc_amp = pick("visc_amp", 50.0)
    p_grad   = pick("p_grad", 4.8)

    x_c     = pick("x_c", 0.012)
    y_0     = pick("y_0", 0.12)
    p_a     = pick("p_a", 1.3)
    p_sigma = pick("p_sigma", 0.3)

    # Anchor (trace-confirmed)
    anchor_logI_acc  = pick("anchor_logI_acc", 0.0250)
    anchor_logI_comp = pick("anchor_logI_comp", 0.3010)

    # --- Structure ---
    xi_z = (1.0 + z_coma) ** (-0.3)

    # --- VSE (engine-faithful) ---
    def eta_acc(x):
        return np.power(np.clip(x / (x + x_c), 0.0, 1.0), p_a)

    def eta_sigma(y):
        return np.power(np.clip(y_0 / (y + y_0), 0.0, 1.0), p_sigma)

    def ratio_r(x, y):
        return x / np.maximum(y + y_0, 1e-30)

    r_ref = float(ratio_r(anchor_logI_acc, anchor_logI_comp))

    r_curr = float(ratio_r(logI_acc, logI_comp))
    eta_grad = float(np.clip((r_curr / max(r_ref, 1e-30)) ** p_grad, 0.0, 1.0))  # cap at 1.0 (no super-eff)

    eta_raw = float(np.clip(eta_acc(logI_acc) * eta_sigma(logI_comp) * eta_grad, 0.0, 2.5))
    eta_ref = float(np.clip(eta_acc(anchor_logI_acc) * eta_sigma(anchor_logI_comp) * 1.0, 1e-30, 2.5))
    eta_eff = float(np.clip(eta_raw / eta_ref, 0.0, 1.0))  # hard guard

    # --- Windows / gates (use engine helpers if present, else mirror known constants) ---
    # acc_gate: center=0.06 width=0.017 (your hard constraint)
    if hasattr(engineV5, "_calculate_acc_gate"):
        acc_gate = float(engineV5._calculate_acc_gate(logI_acc))
    else:
        acc_gate = float(sigmoid((logI_acc - 0.06) / 0.017))

    # W_visc: choose the SAME numbers your engine block uses
    # (example: center=1.2, width=0.15 -> ~0 at 0.6, ~1 at 2.0)
    visc_center = 1.2
    visc_width  = 0.15
    W_visc = float(sigmoid((logI_comp - visc_center) / visc_width))

    if hasattr(engineV5, "_peak_dip_gain"):
        visc_gain = float(engineV5._peak_dip_gain(logI_comp))
    else:
        visc_gain = 1.0

    # Fluid window: if your “direct invariants” path assumes fully fluid-eligible here:
    W_fluid_hi = 1.0

    # --- Enhancements (decoupled assembly) ---
    enh_fluid = float(A0 * W_fluid_hi * eta_eff)
    enh_visc  = float(visc_amp * W_visc * visc_gain * acc_gate)

    enh_total = enh_fluid + enh_visc
    B = 1.0 + enh_total

    v_pred = v_bar * B * xi_z  # IMPORTANT: xi_z applied here, not inside enh

    v_target = 1414.0
    err = abs(v_pred - v_target) / v_target * 100.0

    print(f"Engine constants: A0={A0:.3f}, visc_amp={visc_amp:.2f}, p_grad={p_grad:.2f}, xi_z={xi_z:.4f}")
    print(f"Anchor ratio: r_ref={r_ref:.6e} | Current ratio: r={r_curr:.6e}")
    print(f"VSE: eta_grad={eta_grad:.6f}, eta_eff={eta_eff:.6f}")
    print(f"Visc: W_visc={W_visc:.4f}, visc_gain={visc_gain:.4f}, acc_gate={acc_gate:.4f}")
    print(f"Enh:  enh_fluid={enh_fluid:.4f}, enh_visc={enh_visc:.4f}, enh_total={enh_total:.4f}, B={B:.4f}")
    print("──────────────────────────────────────────────────────────────────────────────────────────")
    print(f"Predicted Velocity: {v_pred:.2f} km/s")
    print(f"Target Velocity:    {v_target:.2f} km/s")
    print(f"Status: {'✅ PASS' if err < 20 else '❌ FAIL'} | Deviation: {err:.1f}%")
    print("══════════════════════════════════════════════════════════════════════════════════════════")

run_isolated_coma_probe()