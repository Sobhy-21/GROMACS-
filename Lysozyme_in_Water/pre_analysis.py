#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os

# Recommended thresholds
TARGET_TEMP = 300.0      # K
TEMP_TOL = 2.0           # K
PRESSURE_TOL = 200.0     # bar, NPT is noisy
POTENTIAL_TOL = 1e4      # arbitrary units

# Define stage properties
STAGE_PROPERTIES = {
    "em": ["Potential"],
    "nvt": ["Temperature", "Potential"],
    "npt": ["Temperature", "Pressure", "Density"],
    "mdrep1": ["Temperature", "Potential"],
    "mdrep2": ["Temperature", "Potential"]
}

def read_xvg(filename):
    """Read XVG file, skip comments, return time and values"""
    data = []
    if not os.path.isfile(filename):
        return None
    with open(filename) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

def analyze_stage(outdir, stage, pdf):
    """Plot and analyze one stage, save to PdfPages"""
    summary = []
    stage_dir = os.path.join(outdir, stage) if stage.startswith("mdrep") else outdir

    if stage.startswith("mdrep"):
        prefix = stage
    else:
        prefix = stage

    for prop in STAGE_PROPERTIES.get(stage, []):
        xvg_file = os.path.join(stage_dir, f"{prefix}_{prop}.xvg")
        data = read_xvg(xvg_file)
        if data is None:
            print(f"Warning: {xvg_file} not found, skipping {prop}.")
            continue

        time = data[:,0]
        values = data[:,1]
        mean = np.mean(values)
        std = np.std(values)

        # Conditional comment
        comment = ""
        ok = True
        if prop == "Temperature":
            if abs(mean - TARGET_TEMP) > TEMP_TOL:
                comment = f"Temperature mean {mean:.2f} K deviates from target {TARGET_TEMP} K."
                ok = False
            else:
                comment = f"Temperature is stable: mean={mean:.2f} ± {std:.2f} K."
        elif prop == "Pressure":
            if std > PRESSURE_TOL:
                comment = f"Pressure fluctuations are large: mean={mean:.2f} ± {std:.2f} bar."
                ok = False
            else:
                comment = f"Pressure is stable: mean={mean:.2f} ± {std:.2f} bar."
        elif prop == "Potential":
            drift = values[-1] - values[0]
            if abs(drift) > POTENTIAL_TOL:
                comment = f"Potential energy shows drift of {drift:.2f} units."
                ok = False
            else:
                comment = f"Potential energy is stable: mean={mean:.2f} ± {std:.2f}."

        # Plot
        plt.figure(figsize=(8,4))
        plt.plot(time, values, label=prop, color='tab:blue')
        plt.xlabel('Time (ps)')
        plt.ylabel(prop)
        plt.title(f'{stage.upper()} - {prop}')
        plt.grid(True)
        plt.tight_layout()

        # Add comment as text box below
        plt.figtext(0.5, -0.1, comment, wrap=True, horizontalalignment='center', fontsize=10)

        # Save to PDF
        pdf.savefig(bbox_inches='tight')
        plt.close()

        # Print
        print(f"{stage} - {prop}: mean={mean:.2f}, std={std:.2f}, ok={ok}")
        print("Comment:", comment)
        summary.append((prop, mean, std, ok, comment))

    return summary

def main():
    # Auto-detect protein outputs in current folder
    output_root = "output"
    for protein in os.listdir(output_root):
        protein_dir = os.path.join(output_root, protein)
        if not os.path.isdir(protein_dir):
            continue

        pdf_file = os.path.join(protein_dir, "pre_analysis.pdf")
        print(f"\n=== Analysis for {protein} → {pdf_file} ===")

        with PdfPages(pdf_file) as pdf:
            # Analyze stages
            stages = ["em", "nvt", "npt", "mdrep1", "mdrep2"]
            for stage in stages:
                analyze_stage(protein_dir, stage, pdf)

        print(f"Saved analysis PDF for {protein}.")

if __name__ == "__main__":
    main()
