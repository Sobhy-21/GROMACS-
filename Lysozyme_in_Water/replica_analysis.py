#!/usr/bin/env python3
"""
Replica MD analysis: RMSD, RMSF, and Radius of Gyration for replica1 and replica2
Outputs a PDF per protein with images + conditional comments
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

# -----------------------------
# User settings
# -----------------------------
output_dir = "output"  # root output directory containing proteins
replicas = ["replica1", "replica2"]

# Stability thresholds (example)
RMSD_tol = 2.0  # Å
Rg_tol = 1.0    # Å, rough check
RMSF_tol = 2.0  # Å

# -----------------------------
# Analysis function
# -----------------------------
def analyze_replica(rep_dir, traj_file, top_file, pdf):
    u = mda.Universe(top_file, traj_file)
    protein = u.select_atoms("protein")

    # Align trajectory
    align.AlignTraj(u, u, select="protein", in_memory=True).run()

    # RMSD
    R = rms.RMSD(u, protein, select="protein").run()
    time = R.rmsd[:,1]
    rmsd_values = R.rmsd[:,2]
    rmsd_mean = np.mean(rmsd_values)
    rmsd_comment = ("RMSD is stable." if rmsd_mean < RMSD_tol 
                    else f"RMSD shows deviations ({rmsd_mean:.2f} Å). Consider longer equilibration or restraints.")

    # RMSF
    protein_rmsf = rms.RMSF(protein).run()
    rmsf_values = protein_rmsf.rmsf
    residue_ids = protein.residues.resids
    rmsf_mean = np.mean(rmsf_values)
    rmsf_comment = ("RMSF is acceptable." if rmsf_mean < RMSF_tol 
                    else f"High RMSF detected ({rmsf_mean:.2f} Å). Check flexible regions or constraints.")

    # Radius of gyration
    rg_values = np.array([protein.radius_of_gyration() for ts in u.trajectory])
    times = np.array([ts.time for ts in u.trajectory])
    rg_mean = np.mean(rg_values)
    rg_comment = ("Radius of gyration stable." if np.std(rg_values) < Rg_tol 
                  else f"Rg fluctuates ({np.std(rg_values):.2f} Å). Might need equilibration adjustment.")

    # Plot to PDF
    def plot_to_pdf(x, y, xlabel, ylabel, title, comment):
        fig, ax = plt.subplots(figsize=(6,4))
        ax.plot(x, y, color="tab:blue")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True)
        # Add comment under figure
        plt.figtext(0.5, -0.05, comment, wrap=True, ha="center", fontsize=9)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close()

    plot_to_pdf(time, rmsd_values, "Time (ps)", "RMSD (Å)", f"{rep_dir} RMSD", rmsd_comment)
    plot_to_pdf(residue_ids, rmsf_values, "Residue ID", "RMSF (Å)", f"{rep_dir} RMSF", rmsf_comment)
    plot_to_pdf(times, rg_values, "Time (ps)", "Radius of Gyration (Å)", f"{rep_dir} Rg", rg_comment)


# -----------------------------
# Main loop
# -----------------------------
def main():
    proteins = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir,d))]
    for protein in proteins:
        protein_dir = os.path.join(output_dir, protein)
        pdf_path = os.path.join(protein_dir, "replicas_analysis.pdf")
        with PdfPages(pdf_path) as pdf:
            for rep in replicas:
                rep_dir = os.path.join(protein_dir, rep)
                if rep == "replica1":
                    traj_file = os.path.join(rep_dir, "mdrep1.xtc")
                    top_file = os.path.join(rep_dir, "mdrep1.tpr")
                elif rep == "replica2":
                    traj_file = os.path.join(rep_dir, "mdrep2.xtc")
                    top_file = os.path.join(rep_dir, "mdrep2.tpr")
                
                if os.path.exists(traj_file) and os.path.exists(top_file):
                    analyze_replica(rep_dir, traj_file, top_file, pdf)
                else:
                    print(f"Skipping {rep_dir}, missing trajectory or topology.")
        print(f"PDF saved: {pdf_path}")


if __name__ == "__main__":
    main()
