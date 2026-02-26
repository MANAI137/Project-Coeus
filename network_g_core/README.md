# ⚙️ Network G Core: The Engine Room V.1.0

This directory contains the immutable mathematical kernel of **Project Coeus**. It serves as the **Active Suspension System** for the metric, calculating the mechanical response of spacetime to baryonic energy density and cosmic expansion.

## 🧱 Core API Components

* **constants.py: The Universal Calibration.** Contains the SI-Grounded thresholds (e.g., 287.35 $M_{\odot}/pc^2$) and unified constants ($\alpha=0.062$, $\gamma=-0.0605$) that govern the global response of the vacuum. This is the "Factory Tune" that maintains the **0.956** $R^2$ performance envelope.

* **engineV5.py: The Active Coilover Logic.** Implements the **Metric Impedance Identity**, providing the mathematical bridge between mass-velocity profiles and gravitational "stiffening." This is the core of the **MagneRide Vacuum** architecture, featuring the **Shatter Wall** (1000 km/s) and **Viscous Gain Shaping** guardrails.

* **voi.py: The Information Mirage.** Manages the Velocity of Information (VoI) calculations used to resolve the **Odometer Paradox** and **Quasar Lensing** anomalies. It handles the "Paving Machine" logic for high-redshift ($z \sim 2$) calculations.



## ⚠️ Kernel Integrity Warning
Modifying `network_g_core/` is strictly prohibited for general validation runs. The equations and constants within this folder represent a **"Correct by Construction"** architecture. They are precisely tuned to maintain the authoritative **18.25 km/s RMSE** across the 175-galaxy SPARC dataset.

> **Even a minor "tweak" to these files can trigger systemic Model Drift, causing the "shocks" to bottom out and rendering the baseline tests in `/tests` invalid.**

## 🛡️ Contribution & Pull Request Policy
This is the "Engine Room." To ensure scientific integrity, we enforce a strict gatekeeping policy for V.1.0:

1.  **Peer-Level Authority:** Pull Requests (PRs) affecting core logic require specialization in General Relativity, Galactic Dynamics, or Computational Physics.
2.  **The Zero-Regression Rule:** Any proposed change must pass the `test_unified.py` suite with **zero regression**. If a change fixes one paradox but breaks the SPARC fit or the **Genzel** recovery, it will be rejected.
3.  **The "Why" Requirement:** Every core modification must be accompanied by a formal technical derivation providing the physical justification for the shift in invariants.



> *"Why is the only real social power; without it, you are powerless. And this is how you come to me—without 'why,' without power."*
> — **The Merovingian**

## 🛠 Usage in Tests
The Core is designed to be a stateless API. When building a new validation module, import the engine directly to calculate the "Active Handling" effect:

```python
from network_g_core import engineV5 as engine

# Calculate the 'Active Suspension' response for a specific mass/velocity profile
# including Pillar 2 (Expansion Damping) and Pillar 4 (Boundary) logic.
v_pred = engine.predict_v_kms(R_kpc, Reff_kpc, Vgas, Vdisk, Vbul, L36, MHI)