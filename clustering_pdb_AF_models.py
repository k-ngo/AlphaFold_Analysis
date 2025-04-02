#!/usr/bin/env python3
"""
clustering_pdb_by_residues_AF_models.py

This script clusters AlphaFold-predicted PDB models by the structural similarity of selected residues.
It computes the all-atom RMSD (using atoms from the specified residues and chains) between each pair of models,
clusters the models using DBSCAN, and outputs cluster details and a bar plot of the average pLDDT prediction confidence score
for each cluster.

For questions, email khoango@ucdavis.edu.
Please cite: https://doi.org/10.7554/eLife.104901.1
"""

import os
import glob
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import PDB
from sklearn.cluster import DBSCAN
from collections import defaultdict
from matplotlib.ticker import AutoMinorLocator
import json
import argparse


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def calculate_all_atom_rmsd_for_residue_ids_in_chain(structure1, structure2, residue_ids, chain_id):
    """
    Calculate the all-atom RMSD between two structures for atoms in the specified residue IDs and chain.
    Only models atoms from the selected chain.
    Returns infinity if the number of atoms does not match.
    """
    super_imposer = PDB.Superimposer()

    atoms1 = [atom for atom in structure1.get_atoms()
              if atom.get_parent().id[1] in residue_ids and atom.get_parent().get_full_id()[2] == chain_id]
    atoms2 = [atom for atom in structure2.get_atoms()
              if atom.get_parent().id[1] in residue_ids and atom.get_parent().get_full_id()[2] == chain_id]

    if len(atoms1) != len(atoms2):
        return float('inf')

    super_imposer.set_atoms(atoms1, atoms2)
    return super_imposer.rms


def cluster_pdb_files(pdb_files, rmsd_threshold, residue_ids, chains):
    """
    Load all PDB structures and compute a pairwise RMSD matrix using the specified residues and chains.
    Cluster the models using DBSCAN (with the given RMSD threshold) and return the clusters as a dict.
    """
    parser = PDB.PDBParser(QUIET=True)
    structures = [parser.get_structure(os.path.basename(pdb_file), pdb_file) for pdb_file in pdb_files]

    pairwise_rmsd = np.zeros((len(structures), len(structures)))
    # Calculate pairwise RMSD by averaging over the specified chains
    for i in range(len(structures)):
        for j in range(i + 1, len(structures)):
            rmsd_values = []
            for chain_id in chains:
                rmsd = calculate_all_atom_rmsd_for_residue_ids_in_chain(structures[i], structures[j], residue_ids, chain_id)
                rmsd_values.append(rmsd)
            avg_rmsd = np.mean(rmsd_values)
            pairwise_rmsd[i, j] = avg_rmsd
            pairwise_rmsd[j, i] = avg_rmsd

    db = DBSCAN(eps=rmsd_threshold, min_samples=1, metric='precomputed')
    db.fit(pairwise_rmsd)

    clusters = defaultdict(list)
    for i, label in enumerate(db.labels_):
        clusters[label].append(pdb_files[i])

    return clusters


def calculate_average_bfactor_for_residue_ids(pdb_file, residue_ids):
    """
    Calculate the average B-factor for atoms in the specified residue IDs.
    In AlphaFold models the B-factor column holds the confidence score.
    Returns None if no atoms are found.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)

    b_factors = [atom.get_bfactor() for atom in structure.get_atoms() if atom.get_parent().id[1] in residue_ids]
    if not b_factors:
        return None
    return np.mean(b_factors)


def save_results_to_json(results, json_file):
    """
    Save clustering results to a JSON file.
    """
    # Convert keys to strings to ensure JSON-serializability.
    results_serializable = {str(key): value for key, value in results.items()}
    with open(json_file, 'w') as file:
        json.dump(results_serializable, file)


def load_results_from_json(json_file):
    """
    Load clustering results from a JSON file.
    """
    with open(json_file, 'r') as file:
        loaded = json.load(file)
    # Convert keys back to int if possible.
    return {int(key): value for key, value in loaded.items()}


def print_cluster_details(clusters):
    """
    Print the details of each cluster, listing the PDB filenames.
    """
    for label, pdb_list in clusters.items():
        print(f"Cluster {label}:")
        pdb_list.sort()
        for pdb_file in pdb_list:
            print(f"- {os.path.basename(pdb_file)}")
        print(' '.join(['"' + os.path.basename(pdb_file) + '"' for pdb_file in pdb_list]))
        print()


def rename_clusters(clusters, residue_ids):
    """
    Rename clusters by sorting them based on the average B-factor (confidence) of the models.
    Cluster labeled 100 is reserved for outliers.
    """
    # Remove outlier cluster (label 100) if present
    outliers = clusters.pop(100, None)
    cluster_labels = list(clusters.keys())
    cluster_avg_b_factors = []
    for cluster_pdb_files in clusters.values():
        avg_vals = []
        for pdb_file in cluster_pdb_files:
            avg = calculate_average_bfactor_for_residue_ids(pdb_file, residue_ids)
            if avg is not None:
                avg_vals.append(avg)
        if avg_vals:
            cluster_avg_b_factors.append(np.mean(avg_vals))
        else:
            cluster_avg_b_factors.append(0)
    # Sort clusters by average B-factor (descending order)
    sorted_clusters = [x for _, x in sorted(zip(cluster_avg_b_factors, cluster_labels), reverse=True)]
    renamed_clusters = {}
    for new_label, old_label in enumerate(sorted_clusters, start=1):
        renamed_clusters[new_label] = clusters[old_label]
    if outliers is not None:
        renamed_clusters[100] = outliers
    return renamed_clusters


def merge_clusters_with_few_members(clusters, min_cluster_size):
    """
    Merge clusters with fewer than min_cluster_size members into a new cluster labeled 100 (outliers).
    """
    new_clusters = {}
    outliers = []
    for label, pdb_list in clusters.items():
        if len(pdb_list) < min_cluster_size:
            outliers.extend(pdb_list)
        else:
            new_clusters[label] = pdb_list
    if outliers:
        new_clusters[100] = outliers
    return new_clusters


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    # ---------------------------
    # ARGUMENT PARSING
    # ---------------------------
    parser = argparse.ArgumentParser(
        description="Cluster PDB models by RMSD over selected residues."
    )
    parser.add_argument("-p", "--pdb_pattern", type=str, required=True,
                        default="alphafold_models/*rank*.pdb",
                        help="(Required) Glob pattern to match predicted PDB files. Example: 'alphafold_models/*rank*.pdb'.")
    parser.add_argument("--residues", type=str, default="227,228,229,230,231",
                        help="Comma-separated list of residue numbers to use for clustering. Example: '227,228,229,230,231'.")
    parser.add_argument("-c", "--chains", type=str, default="A,B,C,D",
                        help="Comma-separated list of chain IDs to process. Example: 'A,B,C,D'.")
    parser.add_argument("--rmsd_threshold", type=float, default=0.35,
                        help="RMSD threshold for clustering (eps value for DBSCAN). Example: 0.35")
    parser.add_argument("--min_cluster_size", type=int, default=3,
                        help="Minimum number of models for a cluster. Clusters with fewer members will be merged into an outlier cluster. Example: 3")
    parser.add_argument("--results_file", type=str, default="pdb_cluster_results.json",
                        help="Filename for saving/loading clustering results (in JSON format).")
    parser.add_argument("--overwrite", action="store_true",
                        help="If specified, overwrite the saved results file.")
    args = parser.parse_args()

    # Print arguments for confirmation.
    print("Arguments:", args)

    # ---------------------------
    # PRE-PROCESSING
    # ---------------------------
    pdb_files = glob.glob(args.pdb_pattern)
    pdb_files = sorted(pdb_files)
    print(f'Number of PDB files: {len(pdb_files)}')
    print('PDB files:', pdb_files)

    # Convert residues and chains arguments to lists.
    residue_ids = [int(x.strip()) for x in args.residues.split(",") if x.strip()]
    chains_to_use = [x.strip() for x in args.chains.split(",") if x.strip()]

    # ---------------------------
    # CLUSTERING
    # ---------------------------
    if (not args.overwrite) and os.path.exists(args.results_file):
        print("Loading clustering results from file...")
        clusters = load_results_from_json(args.results_file)
    else:
        print("Clustering PDB files...")
        clusters = cluster_pdb_files(pdb_files, args.rmsd_threshold, residue_ids, chains_to_use)
        save_results_to_json(clusters, args.results_file)

    clusters = merge_clusters_with_few_members(clusters, args.min_cluster_size)
    renamed_clusters = rename_clusters(clusters, residue_ids)

    # Print cluster details.
    print_cluster_details(renamed_clusters)

    # =============================================================================
    # PLOTTING CLUSTER STATISTICS
    # =============================================================================
    cluster_labels = []
    average_b_factors = []
    std_b_factors = []
    # Remove outlier cluster if empty.
    if 100 in renamed_clusters and not renamed_clusters[100]:
        del renamed_clusters[100]

    for label, pdb_list in renamed_clusters.items():
        print(f"Cluster {label}: ", end='')
        cluster_labels.append("Outliers" if label == 100 else str(label))
        b_vals = [calculate_average_bfactor_for_residue_ids(pdb_file, residue_ids)
                  for pdb_file in pdb_list if calculate_average_bfactor_for_residue_ids(pdb_file, residue_ids) is not None]
        if b_vals:
            avg_b = np.mean(b_vals)
            std_b = np.std(b_vals, ddof=1) if len(b_vals) > 1 else 0
        else:
            avg_b, std_b = 0, 0
        average_b_factors.append(avg_b)
        std_b_factors.append(std_b)
        print(f"{len(pdb_list)} PDB files, Average pLDDT = {round(avg_b, 2)}, Std = {round(std_b, 2)}")

    # Create a bar plot of the clusters using Seaborn.
    sns.set(style='ticks', context='talk')
    plt.figure(figsize=(10, 6))
    palette = sns.color_palette('hls', n_colors=len(cluster_labels) - 1) + [(0.6, 0.6, 0.6)]
    ax = sns.barplot(x=cluster_labels, y=average_b_factors, saturation=1, edgecolor='k', width=0.6, palette=palette)
    plt.xlabel("Clusters")
    plt.ylabel("Average pLDDT")
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
    plt.ylim(0, 100)

    # Annotate bars with average B-factor and number of models.
    for (i, v), std, label in zip(enumerate(average_b_factors), std_b_factors, cluster_labels):
        plt.text(i, v + std + 3, f"{v:.2f}", horizontalalignment='center')
        plt.text(i, v + std + 9.5, f"n={renamed_clusters[int(label) if label.isdigit() else 100] if label == 'Outliers' else len(renamed_clusters[int(label)])}", horizontalalignment='center')
        if len(renamed_clusters[int(label) if label.isdigit() else 100]) > 1:
            plt.errorbar(i, v, yerr=std, color='black', capsize=10, capthick=1, linewidth=1.5)

    sns.despine()
    plt.tight_layout()
    plt.savefig('cluster_pdb.png', dpi=500)
    plt.show()


if __name__ == "__main__":
    main()
