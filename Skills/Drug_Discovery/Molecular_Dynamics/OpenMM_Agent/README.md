# OpenMM Agent (Molecular Dynamics)

**ID:** `biomedical.drug_discovery.molecular_dynamics`
**Version:** 1.0.0
**Status:** Alpha
**Category:** Drug Discovery / Simulation

---

## Overview

The **OpenMM Agent** automates the setup and execution of Molecular Dynamics (MD) simulations. While AlphaFold provides static structures, this agent reveals how proteins *move* and change conformation over time, which is critical for understanding drug binding and allosteric effects.

## Key Capabilities

- **System Setup:** Cleans PDBs, adds hydrogens, solvates in water box, adds ions (neutralization).
- **Simulation Control:** Writes and executes Python scripts for Energy Minimization, Equilibration (NVT/NPT), and Production runs.
- **Analysis:**
    - **RMSD:** Structural stability over time.
    - **RMSF:** Flexibility per residue (identifying binding pockets).

## Integration

Can be chained with **AgentD** or **Docking Agents** to verify ligand stability in the binding pocket.

## References
- *OpenMM Python API*
- *GROMACS*
