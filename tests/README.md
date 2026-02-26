# Project Coeus: Unified Validation Suite V.1.0

This directory serves as the authoritative validation harness for the **Network G Physics Engine**. Each sub-directory represents a fundamental "Pillar of Evidence," testing the engine's consistency across diverse gravitational regimes—from local solar system null tests to high-redshift galactic dynamics.

## 🗂️ Directory Architecture

### [00_unified_test]
**The Master Sync Audit.** A cross-regime sweep verifying that a single set of physical invariants satisfies SPARC, SLACS, and Bullet Cluster dynamics simultaneously without parameter retuning.
* **Key Script:** `test_unified.py`

### [01_sparc_rotation]
**Pillar 2: Galactic Kinematics & Redshift Scaling.**
Validation of the Network G boost across 175 galaxies and early-universe dynamics.
* `test_sparc.py`: Full-sample audit ($R^2 = 0.956$).
* `ngc2841_environmental_influence.py`: Proof of the **10% Neighborhood Stiffening** effect.
* `genzel_paradox.py`: Simulation of high-redshift ($z \sim 2$) rotation curves.
* `plot_diagnostics.py`: Global statistical error analysis and RMSE mapping.



### [02_quasar_lensing]
**Pillar 1: High-Redshift Geometric Lensing.**
Validates vacuum phase enhancement in Strong Lensing systems, aligning the engine with the TDCOSMO/H0LiCOW datasets.
* **Key Script:** `audit_tdcosmo.py`

### [03_atomic_clock]
**Pillar 1: Solar System Shielding (Local Null Test).**
Verification that Network G effects remain below measurable thresholds ($<10^{-10}$) in high-density local environments.
* `test_atomic_clock.py`: Galileo E05 eccentricity validation against GREAT (2018).

### [04_network_gravity]
**Pillar 3: Phase-State & Cluster Dynamics.**
Testing the "Solid Phase" limits and structural discrimination of the engine.
* `test_bullet.py`: Validation of the **Shatter Wall** (1000 km/s boost suppression).
* `test_bullet_ghost.py`: Proof of the **2.33x Boost Ratio** (Baryonic Gas vs. Stellar Structures).
* `test_coma_cluster.py`: Analysis of isolated high-mass clusters.


### [05_galactic_radii]
**Pillar 4: Boundary Logic & Newtonian Baselines.**
Anchoring the engine to the SI-Grounded Phase Threshold.
* `heirarchy.py`: Outside-In resolution anchored at **287.35 $M_{\odot}/pc^2$**.
* `verify_radii_baseline.py`: Pure Newtonian control test (0.1% tolerance).

---

## 🛠️ Global Protocol
All authoritative tests link directly to the core physics engine:
`from network_g_core import engineV5 as engine`

## ⚖️ License & Attribution
All code in this suite is authored by **Miguel Navarro** (ORCID: [0009-0009-5600-7985](https://orcid.org/0009-0009-5600-7985)) and is licensed under **CC BY 4.0**.

---
*Note: Scripts located in `_legacy_archive` are for historical tracking and are not part of the V.1.0 Gospel validation.*