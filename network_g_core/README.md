# ⚙️ Network G Core: The Engine Room V.6 (Hard-Seal)

This directory contains the immutable mathematical kernel of **Project Coeus**. It serves as the **Active Mechanical System** for the metric, calculating the deterministic response of spacetime to baryonic energy density and cosmic expansion.

## 🧱 V.6 Core API Components

* **constants.py: The Universal Calibration.** Contains the SI-Grounded thresholds and unified constants derived from the **Navarro-Fine Structure Identity** ($\alpha_0 \approx 0.06113$). V.6 anchors these values to the **Transverse Dipole Kernel** ($\alpha_{fs} \cdot 8\pi/3$), moving the framework from empirical fit to topological necessity.

* **engineV6.py: The Mechanistic Kernel.** Implements the **Metric Stiffening Identity**, providing the mathematical bridge between mass-energy profiles and gravitational "stiffening." 
    * **Viscous Phase Splice:** Implements the density-dependent $R_{eff}$ threshold ($\rho(R) < \rho_c$).
    * **Shatter Wall EoS:** Features the V.6 stress scalar $\chi \equiv (v^2/c^2) \cdot (\nabla^2 \Phi / \rho_c)$ for automatic GR-relaxation in high-shear environments.



* **voi.py: The Information Mirage.** Manages the Velocity of Information (VoI) calculations. V.6 employs **Gradient Filtering** to ensure the global C-channel impedance (Hubble Tension) does not contaminate the local A-channel potential (Galactic Kinematics).

## ⚠️ Kernel Integrity Warning
Modifying `network_g_core/` is strictly prohibited for general validation runs. The V.6 architecture is **"Correct by Construction."** These equations are locked to maintain the authoritative **18.21 km/s RMSE** across the 175-galaxy SPARC dataset. 

> **Structural Warning:** Even a minor "tweak" to the V.6 kernel can trigger systemic Model Drift, causing the vacuum "shocks" to bottom out and rendering all `/tests` results invalid.

## 🛡️ Contribution & Pull Request Policy
1.  **Peer-Level Authority:** PRs affecting core logic require specialization in General Relativity, Galactic Dynamics, or Computational Physics.
2.  **The Zero-Regression Rule:** Changes must pass the `test_unified.py` suite with **zero regression**. If a change resolves one paradox but breaks the **Genzel Paradox** ($z \approx 2$) recovery, it will be rejected.
3.  **The "Why" Requirement:** Every modification must be accompanied by a formal technical derivation (e.g., a holonomy normalization proof) providing the physical justification for the shift in invariants.



## 🛠 Usage in Tests
The Core is designed as a stateless API:

```python
from network_g_core import engineV6 as engine

# Calculate the V.6 'Active Suspension' response
# Incorporates Pillar 1 (Stiffening), Pillar 2 (Damping), and Pillar 4 (Saturation).
v_pred = engine.predict_v_kms(R_kpc, Reff_kpc, Vgas, Vdisk, Vbul, L36, MHI)