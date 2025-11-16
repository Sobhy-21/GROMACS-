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










