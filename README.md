# GROMACS Automation Scripts

> A collection of automated scripts to streamline GROMACS workflows.

---

## Overview

This repository contains **all the scripts and tools** I use to automate and simplify working with **GROMACS**, the widely used molecular dynamics simulation package.  

Whether you are:
- Preparing input files,  
- Running simulations,  
- Processing trajectories, or  
- Analyzing results,  

these scripts aim to **save time, reduce errors, and standardize workflows**.

---

## Repository Structure

```bash
gromacs/
│
├── prep/                   # Scripts to prepare systems and input files
│   ├── generate_topology.sh
│   └── solvate_system.py
│
├── run/                    # Scripts to automate simulation runs
│   ├── run_md.sh
│   └── batch_submit.py
│
├── analysis/               # Scripts to analyze trajectories and outputs
│   ├── rmsd_analysis.py
│   ├── gyration_plot.sh
│   └── energy_summary.py
│
├── utils/                  # Utility scripts for common tasks
│   ├── file_converter.py
│   └── clean_logs.sh
│
└── README.md               # Repository documentation
