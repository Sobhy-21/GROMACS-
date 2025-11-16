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
│ │ ├── replica2/
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
- Calls `pre_analysis.py` for automated analysis

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
Specifically focuses on replica trajectories:

- Calculates backbone RMSD, C-alpha RMSF, and radius of gyration (Rg)
- Produces plots and PDF summary for further inspection of structural dynamics

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









