# GROMACS Multi-Protein Simulation Pipeline

## Overview
This project provides a fully automated pipeline for preparing, running, and analyzing molecular dynamics (MD) simulations of multiple protein systems in water using GROMACS. The workflow includes energy minimization, NVT and NPT equilibration, replica simulations, and comprehensive post-processing with graphical reports.  

The pipeline is designed to simplify MD simulation tasks, enable reproducibility, and generate high-quality analysis outputs for structural stability assessment.

---

## Directory Structure

```bash
gromacs_tutorials/
│
├── proteins/ # Input protein PDB files
│ ├── protein1.pdb
│ ├── protein2.pdb
│ └── ...
│
├── output/ # Simulation results (auto-created)
│ ├── protein1/
│ │ ├── em.gro
│ │ ├── nvt.gro
│ │ ├── npt.gro
│ │ ├── replica1/
│ │ ├── replica2/ ├── replica_analysis.pdf
│ │ ├── pre_analysis.pdf
│ │ └── plots/
│ └── protein2/
│ └── ...
│
├── multi_protein_pipeline.py # Main script to run full simulations
├── pre_analysis.py # Automated plotting and PDF report generator
└── replica_analysis.py # Replica RMSD/RMSF/Rg analysis

```



---

## Requirements

- [GROMACS](http://www.gromacs.org/) (tested with 2024.2)
- Python 3.12 or higher
- Python libraries:
  - `matplotlib` – for plotting simulation properties
  - `numpy` – for numerical computations
  - `MDAnalysis` – for RMSD/RMSF/Rg analysis
  - `fpdf2` – for PDF report generation

---

## Usage

```bash
module load gromacs
module load python-data/3.12-25.09
python3 multi_protein_pipeline.py proteins
```


## Scripts

### 1. `multi_protein_pipeline.py`
Automates the full simulation workflow:

- Converts PDB files to GROMACS topology (`pdb2gmx`)
- Defines simulation box and solvates the protein
- Adds ions and performs energy minimization (EM)
- Runs NVT and NPT equilibration
- Launches short replica simulations for further sampling
- Calls `pre_analysis.py` and `replica_analysis.py` for automated analysis

### 2. `pre_analysis.py`
Performs post-processing and generates a comprehensive PDF report (`pre_analysis.pdf`) for each protein:

- Reads `.xvg` output files from EM, NVT, NPT, Replica 1, and Replica 2
- Generates high-resolution plots for:
  - Temperature
  - Pressure
  - Density
  - Potential energy
  - RMSD, RMSF, and radius of gyration (replicas)
- Computes averages, fluctuations, and drift
- Evaluates stability and adds conditional comments
- Exports all plots and analysis to a single PDF

### 3. `replica_analysis.py`

Performs detailed analysis of the replica simulations for each protein:

- Reads the trajectory files from Replica 1 and Replica 2
- Calculates structural metrics including:
  - Backbone RMSD
  - C-alpha RMSF
  - Radius of gyration (Rg)
- Generates high-resolution plots for each metric
- Evaluates stability and highlights fluctuations with conditional comments
- Exports all plots and analysis into a single PDF report (`replica_analysis.pdf`) in the protein's output directory


---

## Output

Typical results for a single protein:

| Stage      | Property | Expected Result                              |
|------------|----------|----------------------------------------------|
| EM         | Potential | Stable minimization (low energy)            |
| NVT        | Temperature | ~300 K ± 2 K                              |
| NPT        | Pressure | Stable around 1 bar (±200 bar)               |
| NPT        | Density | Consistent density of water                  |
| Replica 1  | RMSD    | Stable backbone RMSD (~0.1–0.3 nm)           |
| Replica 2  | RMSD    | Stable backbone RMSD (~0.1–0.3 nm)           |
| Replica 1  | RMSF    | Normal flexibility peaks at loops            |
| Replica 2  | RMSF    | Normal flexibility peaks at loops            |
| Replica 1  | Rg      | Stable protein compactness (~1–2 nm)         |
| Replica 2  | Rg      | Stable protein compactness (~1–2 nm)         |



## Customization

Users can modify:

- MDP parameters for EM, NVT, NPT, and replica runs
- Thermostat and barostat settings
- Number of steps, timestep, constraints, or target temperature
- Paths for input/output directories

## Visualization

- Simulation snapshots and structural analyses can be visualized using **PyMOL**, allowing inspection of:
  - Protein conformation
  - Secondary structure
  - Dynamics over time

## Notes

- Each protein is processed independently, enabling parallel execution on HPC clusters.
- All scripts are designed to work on the **Puhti CSC** environment, but can be adapted for other systems with minor adjustments.
- Reports include conditional comments highlighting stability issues and suggestions for further improvement.






