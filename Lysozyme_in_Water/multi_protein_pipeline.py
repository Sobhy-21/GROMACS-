#!/usr/bin/env python3

import os
import sys
import subprocess

GMX = "gmx_mpi"


# ----------------------------
# Utility: run shell commands
# ----------------------------
def run(cmd, input=None):
    print("\n>>>", " ".join(cmd))
    result = subprocess.run(cmd, input=input, text=True)

    if result.returncode != 0:
        print("ERROR running:", " ".join(cmd))
        sys.exit(1)


# ----------------------------
# Write .mdp files
# ----------------------------
def write_mdp_files(outdir):
    mdps = {
        "ions.mdp": """
integrator = steep
emtol = 1000
emstep = 0.01
nsteps = 5000
""",
        "em.mdp": """
integrator  = steep         ; steepest descent
emtol       = 100           ; stop when max force < 100 kJ/mol/nm
emstep      = 0.01          ; step size
nsteps      = 50000         ; maximum number of steps

nstlist     = 1
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
""",
        "nvt.mdp": """
define      = -DPOSRES       ; restrain protein positions
integrator  = md
nsteps      = 50000          ; 100 ps
dt          = 0.002          ; 2 fs
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

constraint_algorithm = lincs
constraints          = h-bonds
lincs_iter           = 1
lincs_order          = 4

cutoff-scheme = Verlet
ns_type      = grid
rcoulomb     = 1.0
rvdw         = 1.0
DispCorr     = EnerPres
coulombtype  = PME
pme_order    = 4
fourierspacing = 0.16

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300

pcoupl      = no
pbc         = xyz

gen_vel     = yes
gen_temp    = 300
gen_seed    = -1

""",
        "npt.mdp": """
define      = -DPOSRES
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

constraint_algorithm = lincs
constraints          = h-bonds
lincs_iter           = 1
lincs_order          = 4

cutoff-scheme = Verlet
ns_type      = grid
rcoulomb     = 1.0
rvdw         = 1.0
DispCorr     = EnerPres
coulombtype  = PME
pme_order    = 4
fourierspacing = 0.16

tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300

pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

pbc         = xyz
gen_vel     = no

""",
        "mdrep1.mdp": """
integrator  = md
nsteps      = 5000      ; 10 ps
dt          = 0.002
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
nstxtcout   = 500
gen-vel     = yes
gen-seed    = 12345
""",
        "mdrep2.mdp": """
integrator  = md
nsteps      = 5000      ; 10 ps
dt          = 0.002
tc-grps     = System
tau-t       = 0.1
ref-t       = 300
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500
nstxtcout   = 500
gen-vel     = yes
gen-seed    = 12346
""",
    }

    for name, text in mdps.items():
        with open(os.path.join(outdir, name), "w") as f:
            f.write(text.strip() + "\n")


# ----------------------------
# extract xvg
# ----------------------------
def extract_xvg(outdir, stage):
    """
    Extract xvg files for analysis.
    Supports:
      - EM: Potential
      - NVT: Temperature, Potential
      - NPT: Temperature, Pressure, Density
      - mdrep1, mdrep2: Temperature, Potential
    """

    # Detect correct .edr file path
    if stage in ["mdrep1", "mdrep2"]:
        edr_file = os.path.join(outdir, stage, f"{stage}.edr")
    else:
        edr_file = os.path.join(outdir, f"{stage}.edr")

    if not os.path.exists(edr_file):
        print(f"Warning: {edr_file} not found, skipping xvg extraction.")
        return

    # Mapping stage → properties
    if stage == "em":
        selection = ["Potential"]

    elif stage == "nvt":
        selection = ["Temperature", "Potential"]

    elif stage == "npt":
        selection = ["Temperature", "Pressure", "Density"]

    elif stage in ["mdrep1", "mdrep2"]:
        selection = ["Temperature", "Potential"]

    else:
        print(f"Unknown stage {stage}, skipping.")
        return

    # Extract each property
    for prop in selection:

        if stage in ["mdrep1", "mdrep2"]:
            xvg_file = os.path.join(outdir, stage, f"{stage}_{prop}.xvg")
        else:
            xvg_file = os.path.join(outdir, f"{stage}_{prop}.xvg")

        print(f"Extracting {prop} → {xvg_file}")

        cmd = [GMX, "energy", "-f", edr_file, "-o", xvg_file]
        run(cmd, input=f"{prop}\n0\n")



# ----------------------------
# Processing a single protein
# ----------------------------
def process_protein(pdb_path, outdir):
    os.makedirs(outdir, exist_ok=True)
    write_mdp_files(outdir)

    # ------------------------------------
    # 1. pdb2gmx
    # ------------------------------------

    # Temporarily switch to the protein's output folder
    pdb_path = os.path.abspath(pdb_path)  # ensure full path
    old_cwd = os.getcwd()
    os.chdir(outdir)

    try:
        run([
            GMX, "pdb2gmx",
            "-ignh",
            "-f", pdb_path,
            "-o", "processed.gro",
            "-p", "topol.top",
            "-ff", "charmm27",
            "-water", "tip3p"
        ], input="0\n")
    finally:
        os.chdir(old_cwd)  # Go back to original working directory

    # ------------------------------------
    # 2. Define box
    # ------------------------------------
    run([
        GMX, "editconf",
        "-f", os.path.join(outdir, "processed.gro"),
        "-o", os.path.join(outdir, "newbox.gro"),
        "-c", "-d", "1.0", "-bt", "cubic"
    ])

    # ------------------------------------
    # 3. Solvate
    # ------------------------------------
    run([
        GMX, "solvate",
        "-cp", os.path.join(outdir, "newbox.gro"),
        "-cs", "spc216.gro",
        "-o", os.path.join(outdir, "solv.gro"),
        "-p", os.path.join(outdir, "topol.top")
    ])


    # ------------------------------------
    # 4. Ions: grompp + genion
    # ------------------------------------
    # Temporarily switch to the protein's output folder
    old_cwd = os.getcwd()
    os.chdir(outdir)
    
    try:
        # grompp for ions
        run([
            GMX, "grompp",
            "-f", "ions.mdp",
            "-c", "solv.gro",
            "-p", "topol.top",
            "-o", "ions.tpr",
            "-maxwarn", "10"
        ])
    
        # genion to add ions
        run([
            GMX, "genion",
            "-s", "ions.tpr",
            "-o", "solv_ions.gro",
            "-p", "topol.top",
            "-neutral"
        ], input="SOL\n")
    
    finally:
        os.chdir(old_cwd)  # return to original working directory



    # ------------------------------------
    # 5. Energy minimization
    # ------------------------------------
    # Temporarily switch to the protein's output folder
    old_cwd = os.getcwd()
    os.chdir(outdir)
    
    try:
        # grompp for EM
        run([
            GMX, "grompp",
            "-f", "em.mdp",
            "-c", "solv_ions.gro",
            "-r", "solv_ions.gro",
            "-p", "topol.top",
            "-o", "em.tpr",
            "-maxwarn", "10"
        ])
    
    
        # mdrun for EM
        run([
            GMX, "mdrun",
            "-deffnm", os.path.join(outdir, "em")
        ])
    
    finally:
        os.chdir(old_cwd)
        
    # extract xvg data for analysis
    extract_xvg(outdir, "em")



    # ------------------------------------
    # 6. NVT
    # ------------------------------------
    # Temporarily switch to the protein's output folder
    old_cwd = os.getcwd()
    os.chdir(outdir)
    
    try:
        # grompp for NVT
        run([
            GMX, "grompp",
            "-f", "nvt.mdp",
            "-c", "em.gro",
            "-r", "em.gro",   # include if using POSRES
            "-p", "topol.top",
            "-o", "nvt.tpr",
            "-maxwarn", "10"
        ])
    
    
        # mdrun for NVT
        run([
            GMX, "mdrun",
            "-deffnm", os.path.join(outdir, "nvt")
        ])
        
    finally:
        os.chdir(old_cwd)
        
        
    # extract xvg data for analysis
    extract_xvg(outdir, "nvt")



    # ------------------------------------
    # 7. NPT
    # ------------------------------------
    # Temporarily switch to the protein's output folder
    old_cwd = os.getcwd()
    os.chdir(outdir)
    
    try:
        # grompp for NPT
        run([
            GMX, "grompp",
            "-f", "npt.mdp",
            "-c", "nvt.gro",
            "-r", "nvt.gro",   # include if using POSRES
            "-p", "topol.top",
            "-o", "npt.tpr",
            "-maxwarn", "10"
        ])
    
    
        # mdrun for NPT
        run([
            GMX, "mdrun",
            "-deffnm", "npt"
        ])
        
    finally:
        os.chdir(old_cwd)
    
    # extract xvg data for analysis
    extract_xvg(outdir, "npt")








    rep1_dir = os.path.join(outdir, "replica1")
    os.makedirs(rep1_dir, exist_ok=True)
    
    # ------------------------------------
    # 8. Replica 1
    # ------------------------------------
    old_cwd = os.getcwd()
    os.chdir(rep1_dir)  # switch to replica1 folder
    
    print("Files before grompp for replica1:", os.listdir(outdir))

    try:
        # grompp for replica 1
        run([
            GMX, "grompp",
            "-f", os.path.join(outdir, "mdrep1.mdp"),  # input MDP in protein folder
            "-c", os.path.join(outdir, "npt.gro"),     # input GRO from protein folder
            "-t", os.path.join(outdir, "npt.cpt"),   # continue from NPT checkpoint
            "-p", os.path.join(outdir, "topol.top"),   # topology in protein folder
            "-o", "mdrep1.tpr",
            "-maxwarn", "10"
        ])
    
        # mdrun for replica 1
        run([
            GMX, "mdrun",
            "-s", "mdrep1.tpr",
            "-deffnm", "mdrep1"  # outputs go here in rep1_dir
        ])
        
        # trjconv: unwrap + center
        run([
            GMX, "trjconv",
            "-s", "mdrep1.tpr",
            "-f", "mdrep1.xtc",
            "-o", "mdrep1_final.xtc",
            "-pbc", "mol",
            "-center"
        ], input="1\n0\n")  # select group 0 for both prompts
        
        
    finally:
        os.chdir(old_cwd)
    
    
    extract_xvg(outdir, "mdrep1")
    

        
        
        


    rep2_dir = os.path.join(outdir, "replica2")
    os.makedirs(rep2_dir, exist_ok=True)
    
    # ------------------------------------
    # 9. Replica 2
    # ------------------------------------
    old_cwd = os.getcwd()
    os.chdir(rep2_dir)  # switch to replica2 folder
    
    
    print("Files before grompp for replica2:", os.listdir(outdir))

    try:
        # grompp for replica 2
        run([
            GMX, "grompp",
            "-f", os.path.join(outdir, "mdrep2.mdp"),  # input MDP in protein folder
            "-c", os.path.join(outdir, "npt.gro"),     # input GRO from protein folder
            "-t", os.path.join(outdir, "npt.cpt"),   # continue from NPT checkpoint
            "-p", os.path.join(outdir, "topol.top"),   # topology in protein folder
            "-o", "mdrep2.tpr",
            "-maxwarn", "10"
        ])
    
        # mdrun for replica 2
        run([
            GMX, "mdrun",
            "-s", "mdrep2.tpr",
            "-deffnm", "mdrep2"  # outputs go here in rep2_dir
        ])
        
        
        # trjconv: unwrap + center
        run([
            GMX, "trjconv",
            "-s", "mdrep2.tpr",
            "-f", "mdrep2.xtc",
            "-o", "mdrep2_final.xtc",
            "-pbc", "mol",
            "-center"
        ], input="1\n0\n")  # select group 0 for both prompts
        
        
        
        
    finally:
        os.chdir(old_cwd)
        
    
    extract_xvg(outdir, "mdrep2")




# ----------------------------
# MAIN
# ----------------------------
"""
def main():
    proteins_dir = "proteins"
    for protein in os.listdir(proteins_dir):
        if protein.endswith(".pdb"):
            pdb_path = os.path.join(proteins_dir, protein)
            outdir = os.path.join("output", protein[:-4])
            print(f"\n===== PROCESSING {protein} =====")
            process_protein(pdb_path, outdir)
"""
def main():
    proteins_dir = "proteins"
    for protein in os.listdir(proteins_dir):
        if protein.endswith(".pdb"):
            pdb_path = os.path.join(proteins_dir, protein)
            outdir = os.path.abspath(os.path.join("output", protein[:-4]))  # use absolute path
            print(f"\n===== PROCESSING {protein} =====")
            print("Output directory:", outdir)
            process_protein(pdb_path, outdir)




if __name__ == "__main__":
    main()




print("Running global MD analysis...")
os.system("python3 pre_analysis.py")


print("Running replica analysis for RMSD/RMSF/Rg...")
os.system("python3 replica_analysis.py")
