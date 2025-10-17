# ---------------- Script Header ----------------
def print_script_header():
    """
    Prints the purpose and description of the ligand-protein MD analysis script.
    """
    header = """
===============================================================================
Molecular Dynamics Ligand-Protein Interaction Analysis Script
===============================================================================

Purpose:
--------
This script performs a comprehensive analysis of a molecular dynamics (MD)
trajectory to evaluate ligand-protein interactions and ligand-solvent 
(water) interactions. It generates plots, computes statistics, and exports
results for residues and waters that interact with the ligand.

Main Features:
--------------
1. Protein Stability:
   - RMSD (Root Mean Square Deviation) of protein backbone
   - RMSF (Root Mean Square Fluctuation) per residue
   - Radius of gyration (Rg)

2. Ligand-Protein Interactions:
   - Residue-ligand contacts (fraction of frames < cutoff, default 5 Å)
   - Hydrogen bonds between ligand and protein residues
   - Other interactions: Hydrophobic, Pi-Pi, Cation-Pi, Salt-Bridge
   - Only residues with interactions above threshold fraction of frames are shown

3. Ligand-Water Interactions:
   - Hydrogen bonds between ligand and water molecules
   - Only waters forming H-bonds in more than a threshold fraction of frames
     (default: 10%) are shown
   - Per-frame counting and classification: Weak, Moderate, Strong

Threshold Rule:
---------------
Only interactions (contacts, H-bonds, hydrophobic, pi-pi, cation-pi, salt-bridges) 
that occur in more than a specified fraction of frames (default: 10%) are shown 
in plots and exported. Very transient interactions below this threshold are ignored.

Outputs:
--------
- Plots for RMSD, RMSF, Rg, residue-ligand contacts, H-bonds, and interaction counts
- CSV file: 'residue_ligand_interactions_filtered.csv' for residues with interactions
- Console summaries with conditional comments for all analyses

Usage:
------
- Set user input parameters: topology_file, trajectory_file, ligand_resname, cutoff, etc.
- Run the script in an environment with MDTraj, NumPy, Matplotlib, and Pandas installed.
- Review plots and printed comments for interpretation of MD results.

Notes:
------
- Distances are measured in nanometers in MDTraj; converted to Å in plots.
- Hydrogen bonds are classified based on distance and angle thresholds.
- Residue labels use chain-resname-resSeq (e.g., CYS-A-126) for clarity.

===============================================================================
"""
    print(header)

# ---------------- Call at the start ----------------
print_script_header()








import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ===== USER INPUT =====
topology_file = "md_rep1_c_final.gro"
trajectory_file = "md_rep1_centered.xtc"
ligand_resname = "UNK"
cutoff = 0.5  # 5 Å in nm
contact_fraction_threshold = 0.1  # consider residue if contact occurs in >10% of frames

# ===== LOAD TRAJECTORY =====
traj = md.load(trajectory_file, top=topology_file)
print("Trajectory loaded:", traj)

# ===== SELECT ATOMS =====
ligand_atoms = traj.topology.select(f"resname {ligand_resname}")
protein_atoms = traj.topology.select("protein")
backbone_atoms = traj.topology.select("protein and backbone")
if len(ligand_atoms) == 0:
    raise ValueError(f"No ligand atoms found for resname '{ligand_resname}'")












# -----------------------------
# 1. RMSD
# -----------------------------
rmsd = md.rmsd(traj, traj, atom_indices=backbone_atoms) * 10  # Å
plt.figure(figsize=(8,5))
plt.plot(traj.time, rmsd)
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (Å)')
plt.title('Protein Backbone RMSD')
plt.grid(True)
plt.show()

avg_rmsd = np.mean(rmsd); std_rmsd = np.std(rmsd)
min_rmsd = np.min(rmsd); max_rmsd = np.max(rmsd)
print(f"RMSD: mean={avg_rmsd:.2f} Å, SD={std_rmsd:.2f} Å, min={min_rmsd:.2f} Å, max={max_rmsd:.2f} Å")
if avg_rmsd < 2.0:
    print("Comment: Protein is stable and equilibrated.")
elif avg_rmsd < 4.0:
    print("Comment: Moderate fluctuations; check flexible loops or termini.")
else:
    print("Comment: High RMSD; protein may be unstable or simulation not equilibrated.")












"""

# -----------------------------
# 2. RMSF
# -----------------------------
rmsf = md.rmsf(traj, traj[0], atom_indices=protein_atoms) * 10  # Å
residues = [atom.residue.index for atom in traj.topology.atoms if atom.index in protein_atoms]
residue_ids = sorted(set(residues))
rmsf_res = [np.mean([rmsf[i] for i, a in enumerate(residues) if a==r]) for r in residue_ids]

plt.figure(figsize=(12,5))
plt.bar(residue_ids, rmsf_res)
plt.xlabel('Residue index'); plt.ylabel('RMSF (Å)')
plt.title('RMSF per Residue')
plt.show()

mean_rmsf = np.mean(rmsf_res); sd_rmsf = np.std(rmsf_res)
high_rmsf_residues = [residue_ids[i] for i, v in enumerate(rmsf_res) if v > 2.5]
print(f"RMSF: mean={mean_rmsf:.2f} Å, SD={sd_rmsf:.2f} Å")
if len(high_rmsf_residues) == 0:
    print("Comment: Protein is rigid and stable.")
else:
    print(f"Comment: High RMSF at residues {high_rmsf_residues} → flexible loops or termini.")


"""

# -----------------------------
# 2. RMSF (with residue name summary)
# -----------------------------
rmsf = md.rmsf(traj, traj[0], atom_indices=protein_atoms) * 10  # Å

# Map each protein atom to its residue
residues = [atom.residue for atom in traj.topology.atoms if atom.index in protein_atoms]
residue_ids = sorted(set([res.index for res in residues]))

# Compute mean RMSF per residue
rmsf_res = []
residue_labels = []
for r in residue_ids:
    res_atoms = [i for i, a in enumerate(residues) if a.index == r]
    rmsf_res.append(np.mean([rmsf[i] for i in res_atoms]))
    res_obj = traj.topology.residue(r)
    residue_labels.append(f"{res_obj.name}{res_obj.resSeq}")

# Plot RMSF (keep numeric x-axis for clarity)
plt.figure(figsize=(12,5))
plt.bar(residue_ids, rmsf_res)
plt.xlabel('Residue index'); plt.ylabel('RMSF (Å)')
plt.title('RMSF per Residue')
plt.show()

# Compute statistics and identify flexible residues
mean_rmsf = np.mean(rmsf_res)
sd_rmsf = np.std(rmsf_res)

# Identify high-flexibility residues (RMSF > 2.5 Å)
high_rmsf_indices = [i for i, v in enumerate(rmsf_res) if v > 2.5]
high_rmsf_residues = [residue_labels[i] for i in high_rmsf_indices]

print(f"RMSF: mean={mean_rmsf:.2f} Å, SD={sd_rmsf:.2f} Å")
if len(high_rmsf_residues) == 0:
    print("Comment: Protein is rigid and stable.")
else:
    print(f"Comment: High RMSF at residues {', '.join(high_rmsf_residues)} → flexible loops or termini.")














# -----------------------------
# 3. Radius of gyration
# -----------------------------
protein_traj = traj.atom_slice(protein_atoms)
rg = md.compute_rg(protein_traj) * 10  # Å
plt.figure(figsize=(8,5))
plt.plot(traj.time, rg)
plt.xlabel('Time (ps)'); plt.ylabel('Radius of Gyration (Å)')
plt.title('Protein Radius of Gyration')
plt.grid(True)
plt.show()

mean_rg = np.mean(rg); sd_rg = np.std(rg); delta_rg = np.max(rg)-np.min(rg)
print(f"Rg: mean={mean_rg:.2f} Å, SD={sd_rg:.2f} Å, delta={delta_rg:.2f} Å")
if delta_rg < 1.0:
    print("Comment: Protein compact and stable.")
else:
    print("Comment: Significant fluctuation → possible unfolding or expansion.")












# -----------------------------
# 4. Residue-ligand contacts (fraction of frames with threshold)
# -----------------------------
threshold_fraction = 0.1  # 10% of frames

# All protein-ligand atom pairs
atom_pairs = np.array([[i, j] for i in ligand_atoms for j in protein_atoms])
distances = md.compute_distances(traj, atom_pairs)  # n_frames x n_pairs

# Map protein atoms to residues
protein_res_indices = [traj.topology.atom(j).residue.index for i,j in atom_pairs]

# Count frames where distance < cutoff
residue_frame_counts = {}
for i, res_idx in enumerate(protein_res_indices):
    contact_frames = np.sum(distances[:, i] < cutoff)
    if res_idx in residue_frame_counts:
        residue_frame_counts[res_idx] += contact_frames
    else:
        residue_frame_counts[res_idx] = contact_frames

# Compute fraction of frames in contact
n_frames = traj.n_frames
residue_contact_fraction = {r: residue_frame_counts[r]/n_frames for r in residue_frame_counts}

# Apply threshold: only keep residues with contact in >= threshold_fraction of frames
residue_contact_fraction = {r: f for r,f in residue_contact_fraction.items() if f >= threshold_fraction}

# Map indices to human-readable labels
residue_labels = {res.index: f"{res.name}-{res.resSeq}" for res in traj.topology.residues if res.index in residue_contact_fraction}

# Classify interactions
interaction_type = {}
for r, f in residue_contact_fraction.items():
    if f < 0.3:
        interaction_type[r] = 'Weak'
    elif f < 0.7:
        interaction_type[r] = 'Moderate'
    else:
        interaction_type[r] = 'Strong'

# Plot only residues above threshold
plt.figure(figsize=(12,5))
colors = {'Weak':'lightblue', 'Moderate':'orange', 'Strong':'red'}
plt.bar([residue_labels[r] for r in residue_contact_fraction.keys()],
        list(residue_contact_fraction.values()),
        color=[colors[interaction_type[r]] for r in residue_contact_fraction])
plt.xlabel('Residue')
plt.ylabel('Fraction of frames in contact')
plt.title(f'Residue-Ligand Contacts (fraction of frames <5 Å, threshold: {threshold_fraction*100:.0f}%)')
plt.xticks(rotation=45, ha='right')
plt.ylim(0, 1)  # fraction from 0 to 1
plt.tight_layout()
plt.show()

# Statistics & summary comment
mean_contact = np.mean(list(residue_contact_fraction.values()))
sd_contact = np.std(list(residue_contact_fraction.values()))
strong_contacts = [residue_labels[r] for r,f in residue_contact_fraction.items() if interaction_type[r]=='Strong']

print(f"Contact frequency: mean={mean_contact:.2f}, SD={sd_contact:.2f}")
if strong_contacts:
    print(f"Comment: Strong contacts at residues {strong_contacts} → likely key ligand-binding residues.")
else:
    print(f"Comment: No strong contacts above {threshold_fraction*100:.0f}% of frames → check ligand placement or simulation time.")




















"""

# --- Other interactions (distance <5 Å) ---
residue_interactions_filtered = {}  # only residues with any interaction
residue_labels_filtered = {}        # store labels like CYS-126

for atom_lig in ligand_atoms:
    for atom_prot in protein_atoms:
        prot_res = traj.topology.atom(atom_prot).residue
        res_idx = prot_res.index
        res_label = f"{prot_res.name}-{prot_res.resSeq}"  # e.g., CYS-126
        dist = np.linalg.norm(traj.xyz[:,atom_lig,:]-traj.xyz[:,atom_prot,:], axis=1)
        if np.any(dist<cutoff):
            if res_idx not in residue_interactions_filtered:
                residue_interactions_filtered[res_idx] = {it:0 for it in interaction_types}
                residue_labels_filtered[res_idx] = res_label
            atom = traj.topology.atom(atom_prot)
            # classify interactions
            if atom.element.symbol=='C' and prot_res.name in ['ALA','VAL','LEU','ILE','MET','PHE','PRO','TRP']:
                residue_interactions_filtered[res_idx]['Hydrophobic'] += 1
            if prot_res.name in ['PHE','TYR','TRP','HIS'] and atom.element.symbol=='C':
                residue_interactions_filtered[res_idx]['Pi-Pi'] += 1
            if prot_res.name in ['LYS','ARG'] and atom.element.symbol=='N':
                residue_interactions_filtered[res_idx]['Cation-Pi'] += 1
            if prot_res.name in ['ASP','GLU','LYS','ARG'] and atom.element.symbol in ['O','N']:
                residue_interactions_filtered[res_idx]['Salt-Bridge'] += 1

# Include H-bonds
for donor_idx, hydrogen_idx, acceptor_idx in hbonds:
    donor_res = traj.topology.atom(donor_idx).residue
    acceptor_res = traj.topology.atom(acceptor_idx).residue
    if donor_res.name==ligand_resname or acceptor_res.name==ligand_resname:
        other_res = acceptor_res if donor_res.name==ligand_resname else donor_res
        res_idx = other_res.index
        res_label = f"{other_res.name}-{other_res.resSeq}"
        if res_idx not in residue_interactions_filtered:
            residue_interactions_filtered[res_idx] = {it:0 for it in interaction_types}
            residue_labels_filtered[res_idx] = res_label
        distances = md.compute_distances(traj, [[donor_idx,acceptor_idx]])
        angles = md.compute_angles(traj, [[donor_idx,hydrogen_idx,acceptor_idx]])*180/np.pi
        strength = classify_hbond(np.mean(distances), np.mean(angles))
        if strength:
            residue_interactions_filtered[res_idx]['H-bond'] += 1

# Prepare plot for residues with any interaction
filtered_residues = list(residue_interactions_filtered.keys())
res_labels = [residue_labels_filtered[r] for r in filtered_residues]
interaction_matrix = {it:[residue_interactions_filtered[r][it] for r in filtered_residues] for it in interaction_types}

# Plot stacked bar chart
x = np.arange(len(filtered_residues))
plt.figure(figsize=(max(10,len(filtered_residues)*0.4),5))
bottom = np.zeros(len(filtered_residues))
colors = {'H-bond':'blue','Hydrophobic':'green','Pi-Pi':'orange','Cation-Pi':'red','Salt-Bridge':'purple'}
for it in interaction_types:
    plt.bar(x, interaction_matrix[it], bottom=bottom, color=colors[it], label=it)
    bottom += np.array(interaction_matrix[it])
plt.xticks(x, res_labels, rotation=90)
plt.xlabel("Residue"); plt.ylabel("Interaction counts")
plt.title("Residue-Ligand Interactions (5 Å) - Only Residues with Any Interaction")
plt.legend()
plt.tight_layout()
plt.show()

# Export CSV
summary_df = pd.DataFrame({
    'Residue': res_labels,
})
for it in interaction_types:
    summary_df[it] = [residue_interactions_filtered[r][it] for r in filtered_residues]
summary_df.to_csv('residue_ligand_interactions_filtered.csv', index=False)
print("Exported residues with any interaction to 'residue_ligand_interactions_filtered.csv'")

"""





# --- Other interactions (distance <5 Å) with threshold ---
residue_interactions_counts = {}  # raw counts per interaction type per residue
residue_labels_filtered = {}      # store labels like CYS-126
n_frames = traj.n_frames
threshold_fraction = 0.1  # 10% of frames
interaction_types = ['H-bond','Hydrophobic','Pi-Pi','Cation-Pi','Salt-Bridge']

# Count interactions per frame
for atom_lig in ligand_atoms:
    for atom_prot in protein_atoms:
        prot_res = traj.topology.atom(atom_prot).residue
        res_idx = prot_res.index
        res_label = f"{prot_res.name}-{prot_res.resSeq}"
        dist = np.linalg.norm(traj.xyz[:, atom_lig, :] - traj.xyz[:, atom_prot, :], axis=1)

        if res_idx not in residue_interactions_counts:
            residue_interactions_counts[res_idx] = {it:0 for it in interaction_types}
            residue_labels_filtered[res_idx] = res_label

        # Count per-frame interactions
        atom = traj.topology.atom(atom_prot)
        if np.any(dist < cutoff):
            if atom.element.symbol=='C' and prot_res.name in ['ALA','VAL','LEU','ILE','MET','PHE','PRO','TRP']:
                residue_interactions_counts[res_idx]['Hydrophobic'] += np.sum(dist < cutoff)
            if prot_res.name in ['PHE','TYR','TRP','HIS'] and atom.element.symbol=='C':
                residue_interactions_counts[res_idx]['Pi-Pi'] += np.sum(dist < cutoff)
            if prot_res.name in ['LYS','ARG'] and atom.element.symbol=='N':
                residue_interactions_counts[res_idx]['Cation-Pi'] += np.sum(dist < cutoff)
            if prot_res.name in ['ASP','GLU','LYS','ARG'] and atom.element.symbol in ['O','N']:
                residue_interactions_counts[res_idx]['Salt-Bridge'] += np.sum(dist < cutoff)

# Include H-bonds
for donor_idx, hydrogen_idx, acceptor_idx in hbonds:
    donor_res = traj.topology.atom(donor_idx).residue
    acceptor_res = traj.topology.atom(acceptor_idx).residue
    if donor_res.name==ligand_resname or acceptor_res.name==ligand_resname:
        other_res = acceptor_res if donor_res.name==ligand_resname else donor_res
        res_idx = other_res.index
        res_label = f"{other_res.name}-{other_res.resSeq}"
        if res_idx not in residue_interactions_counts:
            residue_interactions_counts[res_idx] = {it:0 for it in interaction_types}
            residue_labels_filtered[res_idx] = res_label
        distances = md.compute_distances(traj, [[donor_idx, acceptor_idx]])
        if classify_hbond(np.mean(distances), 180):  # angle approx. 180 for counting
            residue_interactions_counts[res_idx]['H-bond'] += np.sum(distances < cutoff)

# Convert counts → fractions and apply threshold
residues_to_plot = []
interaction_matrix = {it: [] for it in interaction_types}
res_labels_final = []

for res_idx, counts in residue_interactions_counts.items():
    include_res = False
    for it in interaction_types:
        frac = counts[it] / n_frames
        counts[it] = frac  # update counts to fractions
        if frac >= threshold_fraction:
            include_res = True
    if include_res:
        residues_to_plot.append(res_idx)
        res_labels_final.append(residue_labels_filtered[res_idx])
        for it in interaction_types:
            interaction_matrix[it].append(counts[it])

# Plot stacked bar chart with wider layout and space for legend
x = np.arange(len(res_labels_final))
plt.figure(figsize=(max(12,len(res_labels_final)*0.5),6))  # wider figure
bottom = np.zeros(len(res_labels_final))
colors = {'H-bond':'blue','Hydrophobic':'green','Pi-Pi':'orange','Cation-Pi':'red','Salt-Bridge':'purple'}
for it in interaction_types:
    plt.bar(x, interaction_matrix[it], bottom=bottom, color=colors[it], label=it, edgecolor='black', linewidth=0.7)
    bottom += np.array(interaction_matrix[it])

plt.xticks(x, res_labels_final, rotation=90)
plt.xlabel("Residue")
plt.ylabel("Fraction of frames with interaction (0–1)")
plt.ylim(0, 1)  # y-axis from 0 to 1
plt.title(f"Residue-Ligand Interactions (5 Å) - Threshold: {threshold_fraction*100}%")

# Place legend outside the plot on the upper-right corner
plt.legend(title="Interaction type", loc='upper left', bbox_to_anchor=(1,1))

plt.tight_layout(rect=[0,0,0.85,1])  # leave space on the right for legend
plt.show()




















"""

# -----------------------------
# 5. Hydrogen bonds per residue
# -----------------------------
def residue_label(res):
    resseq = getattr(res, 'resSeq', None)
    if resseq is None:
        resseq = getattr(res, 'serial', None)
    chain_id = getattr(res.chain, 'id', None)
    if chain_id is None:
        chain_id = f"{res.chain.index}"
    if resseq is not None:
        return f"{chain_id}-{res.name}-{resseq}"
    else:
        return f"{chain_id}-{res.name}-{res.index}"

# Identify residues close to ligand
close_residues = []
for res in traj.topology.residues:
    if not res.is_protein: continue
    res_atoms = [a.index for a in res.atoms]
    pairs = np.array([[a1, a2] for a1 in res_atoms for a2 in ligand_atoms])
    if pairs.size == 0: continue
    dists = md.compute_distances(traj, pairs)
    if np.min(dists) < cutoff:
        close_residues.append(res)

close_residues = sorted(close_residues, key=lambda r: (getattr(r.chain,'id',r.chain.index), r.index))
close_labels = [residue_label(r) for r in close_residues]

hbonds = md.baker_hubbard(traj, periodic=False)

def classify_hbond(dist_nm, angle_deg):
    dist = dist_nm * 10.0
    if dist <= 2.5 and angle_deg >= 150: return "Strong"
    elif dist <= 3.2 and angle_deg >= 130: return "Moderate"
    elif dist <= 4.0 and angle_deg >= 110: return "Weak"
    else: return None

categories = ["Weak","Moderate","Strong"]
res_counts = {label: {c:0 for c in categories} for label in close_labels}

# Count H-bonds per residue
for donor_idx, hydrogen_idx, acceptor_idx in hbonds:
    donor_res = traj.topology.atom(donor_idx).residue
    acceptor_res = traj.topology.atom(acceptor_idx).residue
    if donor_res.name == ligand_resname or acceptor_res.name == ligand_resname:
        other_res = acceptor_res if donor_res.name == ligand_resname else donor_res
        if other_res in close_residues:
            distances = md.compute_distances(traj, [[donor_idx, acceptor_idx]])
            angles = md.compute_angles(traj, [[donor_idx, hydrogen_idx, acceptor_idx]]) * (180/np.pi)
            strength = classify_hbond(np.mean(distances), np.mean(angles))
            if strength:
                res_counts[residue_label(other_res)][strength] += 1

# Plot H-bonds (all close residues)
weak_counts = [res_counts[label]["Weak"] for label in close_labels]
moderate_counts = [res_counts[label]["Moderate"] for label in close_labels]
strong_counts = [res_counts[label]["Strong"] for label in close_labels]
x = np.arange(len(close_labels))
width = 0.25

plt.figure(figsize=(max(8,len(close_labels)*0.4),5))
plt.bar(x - width, weak_counts, width, label="Weak", color="lightblue")
plt.bar(x, moderate_counts, width, label="Moderate", color="orange")
plt.bar(x + width, strong_counts, width, label="Strong", color="red")
plt.xticks(x, close_labels, rotation=45, ha="right")
plt.xlabel("Residue")
plt.ylabel("H-bond counts")
plt.title(f"Ligand ({ligand_resname}) → Residue H-bonds")
plt.legend()
plt.tight_layout()
plt.show()

# Print summary comment only for residues that actually have H-bonds
res_with_hbonds = {label: counts for label, counts in res_counts.items() if sum(counts.values()) > 0}
total_hbonds = sum([sum(counts.values()) for counts in res_with_hbonds.values()])

print(f"Total H-bonds: {total_hbonds}")
if res_with_hbonds:
    print("\nH-bond summary (only residues with at least one H-bond):")
    for label, counts in res_with_hbonds.items():
        print(f"{label}: Weak={counts['Weak']}, Moderate={counts['Moderate']}, Strong={counts['Strong']}")
else:
    print("Comment: No H-bonds detected → check ligand placement or cutoff")

"""

"""

# -----------------------------
# 5. Hydrogen bonds per residue (with thresholded summary)
# -----------------------------
def residue_label(res):
    resseq = getattr(res, 'resSeq', None)
    if resseq is None:
        resseq = getattr(res, 'serial', None)
    chain_id = getattr(res.chain, 'id', None)
    if chain_id is None:
        chain_id = f"{res.chain.index}"
    if resseq is not None:
        return f"{chain_id}-{res.name}-{resseq}"
    else:
        return f"{chain_id}-{res.name}-{res.index}"

# Identify residues close to ligand
close_residues = []
for res in traj.topology.residues:
    if not res.is_protein: continue
    res_atoms = [a.index for a in res.atoms]
    pairs = np.array([[a1, a2] for a1 in res_atoms for a2 in ligand_atoms])
    if pairs.size == 0: continue
    dists = md.compute_distances(traj, pairs)
    if np.min(dists) < cutoff:
        close_residues.append(res)

close_residues = sorted(close_residues, key=lambda r: (getattr(r.chain,'id',r.chain.index), r.index))
close_labels = [residue_label(r) for r in close_residues]

hbonds = md.baker_hubbard(traj, periodic=False)

def classify_hbond(dist_nm, angle_deg):
    dist = dist_nm * 10.0
    if dist <= 2.5 and angle_deg >= 150: return "Strong"
    elif dist <= 3.2 and angle_deg >= 130: return "Moderate"
    elif dist <= 4.0 and angle_deg >= 110: return "Weak"
    else: return None

categories = ["Weak","Moderate","Strong"]
res_counts = {label: {c:0 for c in categories} for label in close_labels}

# Count H-bonds per residue
for donor_idx, hydrogen_idx, acceptor_idx in hbonds:
    donor_res = traj.topology.atom(donor_idx).residue
    acceptor_res = traj.topology.atom(acceptor_idx).residue
    if donor_res.name == ligand_resname or acceptor_res.name == ligand_resname:
        other_res = acceptor_res if donor_res.name == ligand_resname else donor_res
        if other_res in close_residues:
            distances = md.compute_distances(traj, [[donor_idx, acceptor_idx]])
            angles = md.compute_angles(traj, [[donor_idx, hydrogen_idx, acceptor_idx]]) * (180/np.pi)
            strength = classify_hbond(np.mean(distances), np.mean(angles))
            if strength:
                res_counts[residue_label(other_res)][strength] += 1

# Plot H-bonds (all close residues)
weak_counts = [res_counts[label]["Weak"] for label in close_labels]
moderate_counts = [res_counts[label]["Moderate"] for label in close_labels]
strong_counts = [res_counts[label]["Strong"] for label in close_labels]
x = np.arange(len(close_labels))
width = 0.25

plt.figure(figsize=(max(8,len(close_labels)*0.4),5))
plt.bar(x - width, weak_counts, width, label="Weak", color="lightblue")
plt.bar(x, moderate_counts, width, label="Moderate", color="orange")
plt.bar(x + width, strong_counts, width, label="Strong", color="red")
plt.xticks(x, close_labels, rotation=45, ha="right")
plt.xlabel("Residue")
plt.ylabel("H-bond counts")
plt.title(f"Ligand ({ligand_resname}) → Residue H-bonds")
plt.legend()
plt.tight_layout()
plt.show()

# -----------------------------
# Summary with threshold
# -----------------------------
threshold_fraction = 0.1  # 10% of frames
n_frames = traj.n_frames

res_with_hbonds = {}
for label, counts in res_counts.items():
    total_frames = sum(counts.values()) / n_frames  # fraction of frames
    if total_frames >= threshold_fraction:
        res_with_hbonds[label] = counts

total_hbonds = sum([sum(counts.values()) for counts in res_with_hbonds.values()])

print(f"Total H-bonds (above {threshold_fraction*100}% frames): {total_hbonds}")
if res_with_hbonds:
    print("\nH-bond summary (only residues above threshold):")
    for label, counts in res_with_hbonds.items():
        print(f"{label}: Weak={counts['Weak']}, Moderate={counts['Moderate']}, Strong={counts['Strong']}")
else:
    print(f"Comment: No residues form H-bonds with ligand in ≥{threshold_fraction*100}% of frames")



"""




# -----------------------------
# 5. Hydrogen bonds per residue (with thresholded summary)
# -----------------------------
def residue_label(res):
    resseq = getattr(res, 'resSeq', None)
    if resseq is None:
        resseq = getattr(res, 'serial', None)
    chain_id = getattr(res.chain, 'id', None)
    if chain_id is None:
        chain_id = f"{res.chain.index}"
    if resseq is not None:
        return f"{chain_id}-{res.name}-{resseq}"
    else:
        return f"{chain_id}-{res.name}-{res.index}"

threshold_fraction = 0.1  # 10% threshold
n_frames = traj.n_frames

# Identify residues close to ligand
close_residues = []
for res in traj.topology.residues:
    if not res.is_protein: 
        continue
    res_atoms = [a.index for a in res.atoms]
    pairs = np.array([[a1, a2] for a1 in res_atoms for a2 in ligand_atoms])
    if pairs.size == 0: 
        continue
    dists = md.compute_distances(traj, pairs)
    if np.min(dists) < cutoff:
        close_residues.append(res)

close_residues = sorted(close_residues, key=lambda r: (getattr(r.chain,'id',r.chain.index), r.index))
close_labels = [residue_label(r) for r in close_residues]

# Detect hydrogen bonds
hbonds = md.baker_hubbard(traj, periodic=True)  # enable PBC if applicable

# Helper function to classify strength
def classify_hbond(dist_nm, angle_deg):
    dist = dist_nm * 10.0  # convert nm → Å
    if dist <= 2.5 and angle_deg >= 150: return "Strong"
    elif dist <= 3.2 and angle_deg >= 130: return "Moderate"
    elif dist <= 4.0 and angle_deg >= 110: return "Weak"
    else: return None

categories = ["Weak","Moderate","Strong"]
res_counts = {label: {c:0 for c in categories} for label in close_labels}
res_fractions = {label: 0 for label in close_labels}

# Count H-bonds per residue per frame
for donor_idx, hydrogen_idx, acceptor_idx in hbonds:
    donor_res = traj.topology.atom(donor_idx).residue
    acceptor_res = traj.topology.atom(acceptor_idx).residue
    if donor_res.name == ligand_resname or acceptor_res.name == ligand_resname:
        other_res = acceptor_res if donor_res.name == ligand_resname else donor_res
        if other_res in close_residues:
            distances = md.compute_distances(traj, [[donor_idx, acceptor_idx]])  # shape (n_frames, 1)
            angles = md.compute_angles(traj, [[donor_idx, hydrogen_idx, acceptor_idx]]) * (180/np.pi)
            for frame_idx in range(n_frames):
                strength = classify_hbond(distances[frame_idx,0], angles[frame_idx,0])
                if strength:
                    label = residue_label(other_res)
                    res_counts[label][strength] += 1

# Compute fraction of frames where each residue forms at least one H-bond
res_fractions = {label: sum(counts.values())/n_frames for label, counts in res_counts.items()}

# Filter by threshold (≥10% of frames)
res_above_threshold = {label: counts for label, counts in res_counts.items() if res_fractions[label] >= threshold_fraction}

# Plot only residues above threshold
if res_above_threshold:
    labels = list(res_above_threshold.keys())
    weak_counts = [res_above_threshold[l]["Weak"]/n_frames for l in labels]
    moderate_counts = [res_above_threshold[l]["Moderate"]/n_frames for l in labels]
    strong_counts = [res_above_threshold[l]["Strong"]/n_frames for l in labels]
    x = np.arange(len(labels))
    width = 0.25

    plt.figure(figsize=(max(8,len(labels)*0.4),5))
    plt.bar(x - width, weak_counts, width, label="Weak", color="lightblue")
    plt.bar(x, moderate_counts, width, label="Moderate", color="orange")
    plt.bar(x + width, strong_counts, width, label="Strong", color="red")
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.xlabel("Residue")
    plt.ylabel("Fraction of frames with H-bond (0–1)")
    plt.title(f"Ligand ({ligand_resname}) → Residue H-bonds (≥{threshold_fraction*100:.0f}% frames)")
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.show()
else:
    print(f"No residues form H-bonds with the ligand in ≥{threshold_fraction*100:.0f}% of frames.")

# -----------------------------
# Summary section
# -----------------------------
if res_above_threshold:
    total_hbonds = sum([sum(counts.values()) for counts in res_above_threshold.values()])
    print(f"Total H-bonds (above {threshold_fraction*100:.0f}% frames): {total_hbonds}")
    print("\nResidues forming H-bonds above threshold:")
    for label, counts in res_above_threshold.items():
        print(f"{label}: Weak={counts['Weak']}, Moderate={counts['Moderate']}, Strong={counts['Strong']} (fraction={res_fractions[label]:.2f})")
else:
    print(f"Comment: No residues form H-bonds with ligand in ≥{threshold_fraction*100:.0f}% of frames.")






















# -----------------------------
# Ligand ↔ Water H-bonds (per-frame counting + threshold)
# -----------------------------
contact_fraction_threshold = 0.1  # consider water if H-bond occurs in >10% of frames
categories = ["Weak", "Moderate", "Strong"]
water_counts_raw = {}  # raw counts per frame
n_frames = traj.n_frames

for donor_idx, hydrogen_idx, acceptor_idx in hbonds:
    donor_res = traj.topology.atom(donor_idx).residue
    acceptor_res = traj.topology.atom(acceptor_idx).residue

    # Only ligand <-> water H-bonds
    if donor_res.name == ligand_resname or acceptor_res.name == ligand_resname:
        water_res = acceptor_res if donor_res.name == ligand_resname else donor_res
        if water_res.name.startswith("HOH") or water_res.name.startswith("WAT"):
            # Compute distances and angles per frame
            distances = md.compute_distances(traj, [[donor_idx, acceptor_idx]])  # shape: (n_frames,1)
            angles = md.compute_angles(traj, [[donor_idx, hydrogen_idx, acceptor_idx]]) * (180.0/np.pi)

            label = f"{water_res.name}-{water_res.resSeq}"
            if label not in water_counts_raw:
                water_counts_raw[label] = {c:0 for c in categories}

            # Count per-frame
            for d, a in zip(distances[:,0], angles[:,0]):
                strength = classify_hbond(d, a)
                if strength is not None:
                    water_counts_raw[label][strength] += 1

# Filter waters by contact fraction
water_counts = {}
for label, counts in water_counts_raw.items():
    total_frames_with_hbond = sum(counts.values())
    fraction = total_frames_with_hbond / n_frames
    if fraction >= contact_fraction_threshold:
        water_counts[label] = counts

# ===== PLOT ONLY WATERS ABOVE THRESHOLD =====
if not water_counts:
    print(f"No water molecules form H-bonds with the ligand in more than {contact_fraction_threshold*100}% of frames.")
else:
    water_labels = sorted(water_counts.keys())
    weak_counts = [water_counts[label]["Weak"] for label in water_labels]
    moderate_counts = [water_counts[label]["Moderate"] for label in water_labels]
    strong_counts = [water_counts[label]["Strong"] for label in water_labels]

    x = np.arange(len(water_labels))
    width = 0.25

    plt.figure(figsize=(max(8, len(water_labels)*0.4),5))
    plt.bar(x - width, weak_counts, width, label="Weak")
    plt.bar(x, moderate_counts, width, label="Moderate")
    plt.bar(x + width, strong_counts, width, label="Strong")

    plt.xticks(x, water_labels, rotation=45, ha="right")
    plt.xlabel("Water molecules forming H-bonds")
    plt.ylabel("Number of H-bonds (per-frame count)")
    plt.title(f"Ligand ({ligand_resname}) → Water H-bonds by Strength (>{contact_fraction_threshold*100}% frames)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ===== SUMMARY & CONDITIONAL COMMENTS =====
    print("\nHydrogen bond summary for waters with H-bonds (per-frame counts, above threshold):")
    total_hbonds = 0
    for label in water_labels:
        c = water_counts[label]
        total = c['Weak'] + c['Moderate'] + c['Strong']
        total_hbonds += total
        print(f"{label}: Weak={c['Weak']}, Moderate={c['Moderate']}, Strong={c['Strong']} (Total={total})")
    
    # Conditional comment
    if total_hbonds == 0:
        print("Comment: No significant ligand-water H-bonds observed above threshold.")
    elif total_hbonds <= 5:
        print("Comment: Only a few H-bonds observed → weak solvent interaction.")
    elif total_hbonds <= 15:
        print("Comment: Moderate ligand-water H-bond network → transient interactions.")
    else:
        print("Comment: Strong ligand-water H-bond network → significant solvation around the ligand.")
