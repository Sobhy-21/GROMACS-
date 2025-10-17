
import MDAnalysis as mda
from MDAnalysis.lib import distances
import os
import warnings

# Suppress MDAnalysis PDB warnings
warnings.simplefilter("ignore", UserWarning)

# ================= USER SETTINGS =================
topology = "md_rep1_c_final.gro"       # Your topology
trajectory = "md_rep1_centered.xtc"      # Your trajectory
ligand_sel = "resname UNK"   # Your ligand residue name
residue_ids = [69, 126]      # Lys69, Cys126
dist_cutoff = 3.5             # Å, distance threshold for contacts
out_folder = "hb_frames"     # Folder for snapshots
# =================================================

u = mda.Universe(topology, trajectory)

# Ligand selection
ligand = u.select_atoms(ligand_sel)

# Target residues (CYS126 + LYS69)
protein = u.select_atoms(" or ".join(f"resid {r}" for r in residue_ids))

# Select polar atoms only (N/O)
ligand_polar = ligand.select_atoms("name N* or name O*")
protein_polar = protein.select_atoms("name N* or name O*")

print(f"Ligand polar atoms: {len(ligand_polar)}")
print(f"Target residue polar atoms: {len(protein_polar)}")
print(f"Searching for contacts <= {dist_cutoff:.2f} Å...")

os.makedirs(out_folder, exist_ok=True)

contacts = []  # store detected contacts
frames_with_contacts = set()

for ts in u.trajectory:
    idxs, dists = distances.capped_distance(
        ligand_polar.positions, protein_polar.positions,
        max_cutoff=dist_cutoff, box=u.dimensions
    )

    for (i_l, i_p), d in zip(idxs, dists):
        lig_atom = ligand_polar[i_l]
        prot_atom = protein_polar[i_p]

        # Only keep contacts with CYS126 or LYS69
        if prot_atom.resid not in residue_ids:
            continue

        contacts.append({
            "frame": ts.frame,
            "time_ps": ts.time,
            "lig_atom": lig_atom,
            "prot_atom": prot_atom,
            "distance_A": d
        })
        frames_with_contacts.add(ts.frame)

# Save PDB snapshots for frames with contacts
for f in sorted(frames_with_contacts):
    u.trajectory[f]
    outpdb = os.path.join(out_folder, f"frame_{f:06d}.pdb")
    u.atoms.write(outpdb)

# Print summary
print(f"\nFound {len(contacts)} contacts in {len(frames_with_contacts)} frames.")

if contacts:
    print("\nContacts (first 50 shown):")
    for c in contacts[:50]:
        l = c["lig_atom"]
        p = c["prot_atom"]
        print(f"frame {c['frame']:6d}, time {c['time_ps']:8.2f} ps, "
              f"dist {c['distance_A']:.2f} Å : "
              f"LIG {l.name} ↔ {p.resname}{p.resid}/{p.name}")

print(f"\nPDB snapshots saved in '{out_folder}' for frames with contacts.")
