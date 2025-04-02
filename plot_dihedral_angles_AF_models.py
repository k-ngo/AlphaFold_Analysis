#!/usr/bin/env python3
"""
plot_dihedral_angles_AF_models.py

This script processes AlphaFold predicted PDB models and, optionally, reference PDB models to extract dihedral angles (φ or ψ) for specified residues.
It creates plots of the angle distributions and saves them as PNG files.
For questions, email khoango@ucdavis.edu.
Please cite: https://doi.org/10.7554/eLife.104901.1
"""

import argparse
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import PDB

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_phi_angle(residue):
    """
    Calculate and return the phi dihedral angle (in degrees) for a given residue.
    """
    try:
        c_prev = residue.parent[residue.id[1] - 1]["C"].get_vector()
        n = residue["N"].get_vector()
        ca = residue["CA"].get_vector()
        c = residue["C"].get_vector()
        return PDB.calc_dihedral(c_prev, n, ca, c) * (180.0 / np.pi)
    except Exception:
        return None

def get_psi_angle(residue):
    """
    Calculate and return the psi dihedral angle (in degrees) for a given residue.
    """
    try:
        n = residue["N"].get_vector()
        ca = residue["CA"].get_vector()
        c = residue["C"].get_vector()
        n_next = residue.parent[residue.id[1] + 1]["N"].get_vector()
        return PDB.calc_dihedral(n, ca, c, n_next) * (180.0 / np.pi)
    except Exception:
        return None

def extract_model_number(filename):
    """
    Extract and return the model number from the filename based on the pattern 'rank' followed by digits.
    """
    match = re.search(r'rank[_]?(\d+)', filename, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None

# Mapping from three-letter to one-letter residue codes.
three_to_one = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    # ---------------------------
    # ARGUMENT PARSING
    # ---------------------------
    parser = argparse.ArgumentParser(
        description="Extract dihedral angles from PDB models and plot their distributions."
    )
    # Predicted pattern is required; default example: "alphafold_models/*rank*.pdb"
    parser.add_argument("-p", "--predicted_pattern", type=str, required=True,
                        default="alphafold_models/*rank*.pdb",
                        help="(Required) Glob pattern to match predicted PDB files. Example: 'alphafold_models/*rank*.pdb'.")
    parser.add_argument("-f", "--fig_prefix", type=str,
                        default="dihedral_angles",
                        help="Prefix for output figure file names; each plot will be saved as <prefix>_<angle>.png or combined into one figure if grid layout is provided.")
    parser.add_argument("-r", "--reference_files", type=str,
                        default="o.pdb,i.pdb,ic3.pdb,co.pdb",
                        help="(OPTIONAL) Comma-separated list of reference PDB files. Leave empty (\"\") to skip reference analysis.")
    parser.add_argument("--ref_colors", type=str,
                        default="o.pdb:#5B9BD5,i.pdb:#FC5C24,co.pdb:#7CBB6C,ic3.pdb:blueviolet",
                        help="(OPTIONAL) Comma-separated key:value pairs for reference file colors, e.g., 'o.pdb:#5B9BD5,i.pdb:#FC5C24,...'.")
    parser.add_argument("--ref_marker_styles", type=str,
                        default="o.pdb:o,i.pdb:s,co.pdb:^,ic3.pdb:D",
                        help="(OPTIONAL) Comma-separated key:value pairs for reference file marker styles.")
    parser.add_argument("--ref_file_labels", type=str,
                        default="o.pdb:Open (PDB 5VA2),i.pdb:Inact.-sampling Cluster 2,ic3.pdb:Inact.-sampling Cluster 3,co.pdb:Closed",
                        help="(OPTIONAL) Comma-separated key:value pairs for reference file labels.")
    parser.add_argument("-ds", "--dot_size", type=int, default=100,
                        help="Size of the scatter plot dots.")
    parser.add_argument("-dec", "--dot_edge_color", type=str, default="black",
                        help="Edge color for the scatter plot dots.")
    parser.add_argument("-dpt", "--dot_proximity_threshold", type=float, default=20,
                        help="Minimum angular distance (in degrees) between dots to avoid overlap.")
    parser.add_argument("-c", "--chains", type=str, default="A,B,C,D",
                        help="Comma-separated list of chain IDs to process. Example: 'A,B,C,D'.")
    parser.add_argument("--residues", type=str, default="624,625,626,627,628",
                        help="Comma-separated list of residue numbers to analyze. Example: '624,625,626,627,628'.")
    parser.add_argument("-po", "--predicted_offset", type=int, default=397,
                        help="Offset to apply to residue numbering in predicted models, useful if the reference numbering differs (default models start at 1).")
    parser.add_argument("--angle_type", type=str, choices=["phi", "psi"], default="phi",
                        help="Select which dihedral angle to plot. Choose 'phi' (default) or 'psi'.")
    parser.add_argument("--grid_layout", type=int, nargs=2, metavar=('ROWS', 'COLS'),
                        help="(OPTIONAL) Grid layout for combining all angle plots into one figure. Provide two integers: number of rows and columns (e.g., --grid_layout 2 3).")
    args = parser.parse_args()

    # Print all arguments at startup for confirmation.
    print("Arguments:", args)

    # ---------------------------
    # PRE-PROCESSING ARGUMENTS
    # ---------------------------
    # Retrieve predicted files based on the provided glob pattern.
    predicted_files = glob.glob(args.predicted_pattern)
    # Process reference files if provided.
    reference_files = [x.strip() for x in args.reference_files.split(",") if x.strip()] if args.reference_files else []

    # Process ref_colors into a dictionary.
    ref_colors = {}
    if args.ref_colors:
        for pair in args.ref_colors.split(","):
            key, val = pair.split(":")
            ref_colors[key.strip()] = val.strip()

    # Process ref_marker_styles into a dictionary.
    ref_marker_styles = {}
    if args.ref_marker_styles:
        for pair in args.ref_marker_styles.split(","):
            key, val = pair.split(":")
            ref_marker_styles[key.strip()] = val.strip()

    # Process ref_file_labels into a dictionary.
    ref_file_labels = {}
    if args.ref_file_labels:
        for pair in args.ref_file_labels.split(","):
            key, val = pair.split(":", 1)
            ref_file_labels[key.strip()] = val.strip()

    chains_to_use = [x.strip() for x in args.chains.split(",") if x.strip()]
    # Convert the comma-separated list of residues to a list of integers.
    residues = [int(x.strip()) for x in args.residues.split(",") if x.strip()]
    predicted_offset = args.predicted_offset
    angle_type = args.angle_type.lower()
    # Use capitalized angle type for keys (e.g., "Phi" or "Psi").
    angle_prefix = "Phi" if angle_type == "phi" else "Psi"

    # =============================================================================
    # PROCESS PREDICTED MODELS
    # =============================================================================
    # "pred" is an arbitrary identifier used to label the parsed predicted structure.
    pdb_parser = PDB.PDBParser(QUIET=True)
    predicted_rows = []
    for pdb_file in predicted_files:
        structure = pdb_parser.get_structure("pred", pdb_file)
        model_num = extract_model_number(pdb_file)
        for model in structure:
            for chain in model:
                if chain.id not in chains_to_use:
                    continue
                # Initialize a row for each valid chain.
                row = {"PDB_File": pdb_file, "ModelNumber": model_num, "Chain": chain.id}
                valid = True
                # Process each specified residue for the chosen angle type.
                for r in residues:
                    try:
                        # Adjust residue number for predicted models.
                        res = chain[(" ", r + predicted_offset, " ")]
                        if angle_type == "phi":
                            row[f"Phi_{r}"] = get_phi_angle(res)
                        else:
                            row[f"Psi_{r}"] = get_psi_angle(res)
                    except Exception:
                        valid = False
                if valid:
                    predicted_rows.append(row)
    df_pred = pd.DataFrame(predicted_rows)

    # =============================================================================
    # PROCESS REFERENCE MODELS (OPTIONAL)
    # =============================================================================
    if reference_files:
        reference_rows = []
        # "ref" is an arbitrary identifier used to label the parsed reference structure.
        for pdb_file in reference_files:
            structure = pdb_parser.get_structure("ref", pdb_file)
            # Use user-provided label; if not available, fallback to the filename.
            label = ref_file_labels.get(pdb_file, pdb_file)
            for model in structure:
                for chain in model:
                    if chain.id not in chains_to_use:
                        continue
                    row = {"PDB_File": pdb_file, "Chain": chain.id, "State": label}
                    valid = True
                    for r in residues:
                        try:
                            # For reference models, no offset is applied.
                            res = chain[(" ", r, " ")]
                            if angle_type == "phi":
                                row[f"Phi_{r}"] = get_phi_angle(res)
                            else:
                                row[f"Psi_{r}"] = get_psi_angle(res)
                        except Exception:
                            valid = False
                    if valid:
                        reference_rows.append(row)
        df_ref = pd.DataFrame(reference_rows)
    else:
        df_ref = pd.DataFrame()

    # =============================================================================
    # GET RESIDUE LABELS
    # =============================================================================
    # For each residue, obtain a label (one-letter code with residue number).
    # If reference models are provided, use the first reference model; otherwise, use the first predicted model.
    residue_labels = {}
    if reference_files:
        # Use the first available reference model.
        for pdb_file in reference_files:
            structure = pdb_parser.get_structure("ref", pdb_file)
            for model in structure:
                for chain in model:
                    if chain.id in chains_to_use:
                        for r in residues:
                            try:
                                res = chain[(" ", r, " ")]
                                resname = res.get_resname()
                                residue_labels[r] = f"{three_to_one.get(resname.upper(), '?')}{r}"
                            except Exception:
                                residue_labels[r] = f"R{r}"
                        break
                break
            break
    else:
        # No reference provided; use the first predicted model.
        if predicted_files:
            structure = pdb_parser.get_structure("pred", predicted_files[0])
            for model in structure:
                for chain in model:
                    if chain.id in chains_to_use:
                        for r in residues:
                            try:
                                res = chain[(" ", r - predicted_offset, " ")]
                                resname = res.get_resname()
                                residue_labels[r] = f"{three_to_one.get(resname.upper(), '?')}{r}"
                            except Exception:
                                residue_labels[r] = f"R{r}"
                        break
                break
        else:
            for r in residues:
                residue_labels[r] = f"R{r}"

    # =============================================================================
    # PLOTTING
    # =============================================================================
    # Create a list of angle columns, e.g., ["Phi_624", "Phi_625", ...] or ["Psi_624", "Psi_625", ...]
    angle_columns = [f"{angle_prefix}_{r}" for r in residues]

    if args.grid_layout:
        # Combined figure: use specified grid layout.
        nrows, ncols = args.grid_layout
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4.5 * ncols, 3.5 * nrows))
        axes = axes.flatten()
        for i, angle in enumerate(angle_columns):
            ax = axes[i]
            # Combine predicted and reference values (if available) for the current angle.
            if not df_ref.empty:
                values = pd.concat([df_pred, df_ref])[angle].dropna()
            else:
                values = df_pred[angle].dropna()
            bins = np.arange(-180, 190, 10)
            hist_vals, bin_edges, _ = ax.hist(values, bins=bins, color="lightgray", edgecolor="black", alpha=0.85, label="AF models")
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            # Determine label text based on residue.
            rnum = angle.split("_")[1]
            label_text = f"{residue_labels.get(int(rnum), 'R'+rnum)} {angle_prefix}"
            ax.set_xlabel(f"{label_text} Angle (degrees)", fontsize=12)
            ax.set_ylabel("Count", fontsize=12)
            ax.set_title(f"Residue {rnum}, {angle_prefix} Distribution")

            # Plot reference dots if reference data exists.
            if not df_ref.empty:
                placed_dots = []
                for ref_file in reference_files:
                    ref_vals = df_ref[df_ref["PDB_File"] == ref_file][angle].dropna()
                    color = ref_colors.get(ref_file, 'black')
                    marker = ref_marker_styles.get(ref_file, 'o')
                    state_label = ref_file_labels.get(ref_file, ref_file)
                    for val in ref_vals:
                        closest_idx = np.argmin(np.abs(bin_centers - val))
                        center_val = bin_centers[closest_idx]
                        base_y = hist_vals[closest_idx]
                        if base_y == 0:
                            base_y = np.max(hist_vals) * 0.05
                        y_offset = 0
                        while any(abs(center_val - x0) <= args.dot_proximity_threshold and abs(base_y + y_offset - y0) < args.dot_size / 10 for x0, y0 in placed_dots):
                            y_offset += args.dot_size / 10
                        ax.scatter(center_val, base_y + args.dot_size / 5 + y_offset, color=color, marker=marker,
                                   s=args.dot_size, edgecolor=args.dot_edge_color, linewidth=1.2, alpha=0.8, label=state_label)
                        placed_dots.append((center_val, base_y + y_offset))
                handles, labels = ax.get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                ax.legend(by_label.values(), by_label.keys(), fontsize=9, frameon=False)
        # Turn off any unused subplots.
        for j in range(len(angle_columns), len(axes)):
            axes[j].axis('off')
        plt.tight_layout()
        output_filename = f"{args.fig_prefix}_combined.png"
        plt.savefig(output_filename, dpi=400)
        plt.close(fig)
    else:
        # Plot each angle in a separate figure.
        for angle in angle_columns:
            fig, ax = plt.subplots(figsize=(6, 4))
            if not df_ref.empty:
                values = pd.concat([df_pred, df_ref])[angle].dropna()
            else:
                values = df_pred[angle].dropna()
            bins = np.arange(-180, 190, 10)
            hist_vals, bin_edges, _ = ax.hist(values, bins=bins, color="lightgray", edgecolor="black", alpha=0.85, label="AF models")
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            rnum = angle.split("_")[1]
            label_text = f"{residue_labels.get(int(rnum), 'R'+rnum)} {angle_prefix}"
            ax.set_xlabel(f"{label_text} Angle (degrees)", fontsize=12)
            ax.set_ylabel("Count", fontsize=12)
            ax.set_title(f"Residue {rnum}, {angle_prefix} Distribution")

            if not df_ref.empty:
                placed_dots = []
                for ref_file in reference_files:
                    ref_vals = df_ref[df_ref["PDB_File"] == ref_file][angle].dropna()
                    color = ref_colors.get(ref_file, 'black')
                    marker = ref_marker_styles.get(ref_file, 'o')
                    state_label = ref_file_labels.get(ref_file, ref_file)
                    for val in ref_vals:
                        closest_idx = np.argmin(np.abs(bin_centers - val))
                        center_val = bin_centers[closest_idx]
                        base_y = hist_vals[closest_idx]
                        if base_y == 0:
                            base_y = np.max(hist_vals) * 0.05
                        y_offset = 0
                        while any(abs(center_val - x0) <= args.dot_proximity_threshold and abs(base_y + y_offset - y0) < args.dot_size / 10 for x0, y0 in placed_dots):
                            y_offset += args.dot_size / 10
                        ax.scatter(center_val, base_y + args.dot_size / 5 + y_offset, color=color, marker=marker,
                                   s=args.dot_size, edgecolor=args.dot_edge_color, linewidth=1.2, alpha=0.8, label=state_label)
                        placed_dots.append((center_val, base_y + y_offset))
                handles, labels = ax.get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                ax.legend(by_label.values(), by_label.keys(), fontsize=9, frameon=False)

            plt.tight_layout()
            output_filename = f"{args.fig_prefix}_{angle}.png"
            plt.savefig(output_filename, dpi=400)
            plt.close(fig)

if __name__ == "__main__":
    main()
