# ⚙️ Project Coeus: Unified Validation Suite (V.1.0)

This directory serves as the authoritative validation harness for the **Network G Physics Engine**. Each sub-directory represents a fundamental "Pillar of Evidence," verifying the engine's consistency across diverse gravitational regimes—from local solar system null tests to the high-redshift "Impossible Galaxy" assembly.

## 🗂️ Directory Architecture

### [00_unified_test]
**The Master Sync Audit.** A cross-regime sweep verifying that a single set of physical invariants satisfies SPARC, SLACS, and Bullet Cluster dynamics simultaneously without per-galaxy parameter retuning.
* **Key Script:** `test_unified.py`

### [01_sparc_rotation]
**Pillar 2: Galactic Kinematics & Redshift Scaling.**
Validation of the Network G boost across the 175-galaxy SPARC sample and early-universe dynamics.
* `test_sparc.py`: Full-sample audit ($R^2 = 0.9563$).
* `ngc2841_environmental_influence.py`: Proof of the **Network Neighbor** influence in high-mass spirals.
* `genzel_paradox.py`: Simulation of high-redshift ($z \sim 2$) rotation curves using expansion damping logic.
* `plot_diagnostics.py`: Global statistical error analysis and RMSE mapping ($\pm 18.21$ km/s).

### [02_quasar_lensing]
**Pillar 1: High-Redshift Geometric Lensing.**
Validates vacuum phase enhancement in Strong Lensing systems, aligning the engine with the TDCOSMO/H0LiCOW datasets.
* **Key Script:** `audit_tdcosmo.py`

### [03_atomic_clock]
**Pillar 1: Solar System Shielding (Local Null Test).**
Verification of the screening threshold, ensuring Network G effects remain below measurable limits in high-density local environments.
* `test_atomic_clock.py`: Galileo E05 eccentricity validation against GREAT (2018).

### [04_network_gravity]
**Pillar 3: Phase-State & Cluster Dynamics.**
Testing the mechanical limits of the vacuum lattice and structural phase discrimination.
* `test_bullet.py`: Validation of the **Shatter Wall** (boost suppression in high-velocity collisions).
* `test_bullet_ghost.py`: Proof of the **2.33x Boost Ratio** (Baryonic Gas vs. Stellar Structures).
* `test_coma_cluster.py`: Analysis of the transition in high-mass clusters.

### [05_galactic_radii]
**Pillar 4: Boundary Logic & Newtonian Baselines.**
Anchoring the engine to the SI-Grounded Phase Threshold and the Outside-In resolution protocol.
* `hierarchy.py`: Analysis of mass-density hierarchies.
* `verify_radii_baseline.py`: Pure Newtonian control test (0.1% tolerance).

### [06_law_vi]
**Pillar 1 & 2: The Law of Vacuum Impedance.**
Validation of the photon-sector distance law and the resolution of the Hubble Tension.
* **Key Script:** `test_hubble_normalization.py`

---

## 🛠️ Global Protocol
All authoritative tests link directly to the core physics engine:
`from network_g_core import engineV6 as engine`

## ⚖️ License & Attribution
All code in this suite is authored by **Miguel Antonio Navarro** (ORCID: [0009-0009-5600-7985](https://orcid.org/0009-0009-5600-7985)) and is licensed under **CC BY 4.0**.