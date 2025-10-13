#!/bin/bash
# =========================================================
# GROMACS property extraction script
# Extracts Temperature, Pressure, and Potential Energy
# from an .edr file and saves as .xvg files
# =========================================================

# Usage:
#   bash extract_props.sh md_rep1.edr
# =========================================================

# Input energy file (default = md.edr)
EDR_FILE=${1:-md_rep1.edr}

# Check if file exists
if [ ! -f "$EDR_FILE" ]; then
  echo "Error: Energy file '$EDR_FILE' not found!"
  echo "Usage: bash extract_props.sh md_rep1.edr"
  exit 1
fi

echo "Extracting from $EDR_FILE ..."

# Temperature
echo "Temperature" | gmx_mpi energy -f "$EDR_FILE" -o run_temp.xvg
# Pressure
echo "Pressure" | gmx_mpi energy -f "$EDR_FILE" -o run_press.xvg
# Potential energy
echo "Potential" | gmx_mpi energy -f "$EDR_FILE" -o run_poten.xvg

echo "Extraction complete."
echo "Generated files:"
ls -1 run_*.xvg
