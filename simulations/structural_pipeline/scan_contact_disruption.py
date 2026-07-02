"""
scan_contact_disruption.py — Ranks candidate second-shell mutations
in the E. coli NemA FMN pocket (PDB 8BPQ, chain A) by a dimensionless
contact-disruption score.

This is NOT a kcal/mol Delta-Delta-G calculation. It does not use
FoldX, Rosetta, or any thermodynamic energy function. It is a
qualitative screening ranker based on the fraction of heavy-atom
contacts that would be lost if a candidate residue were
in-silico truncated to alanine (i.e., the side chain replaced by a
methyl group at C-beta). The methodological basis is the contact-
cutoff analysis used in protein engineering (Kannan & Vishveshwara
1999, J Mol Biol 292(2):441-464; verified via PubMed PMID 10493887)
and in the Rosetta community's hotspot-residue screening.

For each candidate residue X:
  1. C_X = number of protein heavy atoms within 5.0 A of any X
          heavy atom (excluding X's own atoms).
  2. C_X_ala = same, but with X's side chain truncated to C-beta
              (the alanine-like form).
  3. Score = (C_X - C_X_ala) / C_X
  Score in [0, 1]:
    - < 0.3: residue is structurally tolerant of -> Ala
             (few contacts lost, mutation likely accommodated)
    - 0.3 to 0.6: intermediate
    - > 0.6: residue is structurally costly to mutate
             (many contacts lost, mutation likely destabilizing)

A weighted variant uses 1/distance weights to favor tight contacts.

Inputs:
  - simulations/structural_pipeline/data/8BPQ_chainA.pdb
  - simulations/results/nema_active_site.json (from analyze_active_site.py)

Outputs:
  - simulations/results/nema_contact_disruption.csv
  - simulations/results/nema_contact_disruption.png
"""

import csv
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import iter_heavy_atoms, iter_residue_heavy_coords, side_chain_trim_to_ala

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fetch_nemA_structure import PDB_ID, fetch_8bpq
from analyze_active_site import (
    CHAIN_A_PDB, FMN_POCKET_CUTOFF, get_chain_a_residues,
    pocket_residues_from_fmn, AROMATIC_HYDROPHOBIC,
)


# Tunable contact-cutoff for "is this contact contributing to the
# residue's local environment"
CONTACT_CUTOFF_A = 5.0

# Tolerance bands for the score interpretation
SCORE_TOLERANT = 0.3
SCORE_COSTLY = 0.6


def residue_heavy_coords(residue):
    """Return (atom_names, coords) for heavy atoms of a residue.

    For the in-silico Ala variant (side chain trimmed to CB), only
    backbone (N, CA, C, O) plus CB are kept.
    """
    keep = side_chain_trim_to_ala(residue)
    names = []
    coords = []
    for atom in iter_heavy_atoms(residue):
        if atom.name.strip() in keep:
            names.append(atom.name.strip())
            coords.append(atom.get_vector().get_array())
    return names, np.array(coords)


def full_heavy_coords(residue):
    """Return coords of all heavy atoms in a residue (for the WT form)."""
    coords = []
    for atom in iter_heavy_atoms(residue):
        coords.append(atom.get_vector().get_array())
    return np.array(coords)


def count_contacts(coord_set_a, coord_set_b, cutoff=CONTACT_CUTOFF_A,
                   exclude_pairs=None):
    """Count pairs (a_i, b_j) with distance < cutoff.

    exclude_pairs: optional set of frozenset({a_idx, b_idx}) to skip
    (used to exclude X's own atoms from the count).
    """
    if exclude_pairs is None:
        exclude_pairs = set()
    if len(coord_set_a) == 0 or len(coord_set_b) == 0:
        return 0, 0.0
    a = np.array(coord_set_a)
    b = np.array(coord_set_b)
    diff = a[:, None, :] - b[None, :, :]
    d = np.sqrt(np.sum(diff * diff, axis=2))
    n_contacts = int((d < cutoff).sum())
    # Weighted sum: sum of (1/d) over contacts
    contact_mask = d < cutoff
    if contact_mask.any():
        w = np.where(contact_mask, 1.0 / np.maximum(d, 0.1), 0.0)
        weighted = float(w.sum())
    else:
        weighted = 0.0
    return n_contacts, weighted


def compute_contact_disruption_score(residue, protein_residues, cutoff=CONTACT_CUTOFF_A):
    """Compute the contact-disruption score for a single residue.

    Returns (score, weighted_score, n_contacts_wt, n_contacts_ala).
    """
    # 1. WT heavy-atom coordinates of this residue
    wt_coords = full_heavy_coords(residue)

    # 2. In-silico Ala coordinates (backbone + CB only)
    ala_coords = residue_heavy_coords(residue)[1]

    # 3. For every OTHER residue, compute the heavy-atom coords and
    #    count contacts to the WT side and to the Ala side
    n_wt = 0
    n_ala = 0
    w_wt = 0.0
    w_ala = 0.0
    for other in protein_residues:
        if other is residue:
            continue
        other_coords = full_heavy_coords(other)
        if len(other_coords) == 0:
            continue
        n_w, w_w = count_contacts(wt_coords, other_coords, cutoff=cutoff)
        n_a, w_a = count_contacts(ala_coords, other_coords, cutoff=cutoff)
        n_wt += n_w
        n_ala += n_a
        w_wt += w_w
        w_ala += w_a

    if n_wt == 0:
        return None, None, 0, 0
    score = (n_wt - n_ala) / n_wt
    weighted = (w_wt - w_ala) / w_wt if w_wt > 0 else None
    return float(score), (float(weighted) if weighted is not None else None), \
           int(n_wt), int(n_ala)


def scan_candidates(structure_path, candidates, all_residues):
    """Compute the contact-disruption score for each candidate residue.

    Returns a list of dicts with: residue_id, residue_name, chain,
    n_contacts_wt, n_contacts_ala, contact_disruption_score,
    weighted_score, score_band.
    """
    out = []
    for c in candidates:
        score, weighted, n_wt, n_ala = compute_contact_disruption_score(
            c, all_residues, cutoff=CONTACT_CUTOFF_A
        )
        if score is None:
            continue
        # Classify the score into a tolerance band
        if score < SCORE_TOLERANT:
            band = "tolerant"
        elif score > SCORE_COSTLY:
            band = "costly"
        else:
            band = "intermediate"
        out.append({
            "residue_id": c.id[1],
            "residue_name": c.get_resname().strip(),
            "chain": c.get_parent().id,
            "n_contacts_wt": n_wt,
            "n_contacts_ala": n_ala,
            "contact_disruption_score": round(score, 4),
            "weighted_score": (round(weighted, 4) if weighted is not None else ""),
            "score_band": band,
        })
    out.sort(key=lambda x: x["contact_disruption_score"])
    return out


def write_csv(results, out_path):
    """Write the scan results to CSV."""
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def plot_scores(results, out_path):
    """2-panel figure:
      (a) Bar chart of contact-disruption scores, sorted ascending.
          Bars colored by tolerance band.
      (b) Weighted vs unweighted score scatter, labeled by residue.
    """
    sorted_r = sorted(results, key=lambda x: x["contact_disruption_score"])
    labels = [f"{r['residue_name']}{r['residue_id']}" for r in sorted_r]
    scores = [r["contact_disruption_score"] for r in sorted_r]
    weighted = [r["weighted_score"] for r in sorted_r if r["weighted_score"] != ""]
    band_colors = {
        "tolerant": "#2ca02c",
        "intermediate": "#ff7f0e",
        "costly": "#d62728",
    }
    colors = [band_colors[r["score_band"]] for r in sorted_r]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel (a)
    axes[0].barh(labels, scores, color=colors, edgecolor="black", linewidth=0.5)
    axes[0].axvline(SCORE_TOLERANT, color="#2ca02c", linestyle="--",
                    lw=1.5, label=f"tolerant cutoff ({SCORE_TOLERANT})")
    axes[0].axvline(SCORE_COSTLY, color="#d62728", linestyle="--",
                    lw=1.5, label=f"costly cutoff ({SCORE_COSTLY})")
    axes[0].set_xlabel("contact-disruption score (dimensionless)")
    axes[0].set_title("(a) Contact-disruption scores for second-shell candidates\n"
                      "(lower = more tolerant of -> Ala substitution)")
    axes[0].invert_yaxis()
    axes[0].grid(True, axis="x", alpha=0.3)
    axes[0].legend(loc="lower right", fontsize=9)

    # Panel (b)
    if weighted and len(weighted) == len(scores):
        axes[1].scatter(scores, weighted, s=80, c=colors, edgecolor="black",
                        linewidth=0.7, zorder=3)
        for r, x, y in zip(sorted_r, scores, weighted):
            axes[1].annotate(
                f"{r['residue_name']}{r['residue_id']}",
                (x, y), xytext=(4, 4), textcoords="offset points",
                fontsize=8,
            )
        axes[1].axvline(SCORE_TOLERANT, color="#2ca02c", linestyle="--", lw=1.5)
        axes[1].axvline(SCORE_COSTLY, color="#d62728", linestyle="--", lw=1.5)
        axes[1].set_xlabel("contact-disruption score (count-based)")
        axes[1].set_ylabel("weighted score (1/distance weights)")
        axes[1].set_title("(b) Count-based vs distance-weighted scores")
        axes[1].grid(True, alpha=0.3)
    else:
        axes[1].text(0.5, 0.5, "No weighted scores available",
                     ha="center", va="center", transform=axes[1].transAxes)

    fig.suptitle(
        f"Contact-disruption scan of NemA second-shell candidates (PDB {PDB_ID})",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    # Load structure
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    cif, pdb = fetch_8bpq(data_dir)
    if pdb is None or not os.path.exists(pdb):
        # Extract from cif as before
        from Bio.PDB import MMCIFParser, PDBIO, Select
        parser = MMCIFParser(QUIET=True)
        s = parser.get_structure(PDB_ID, cif)

        class ChainASelect(Select):
            def accept_chain(self, chain):
                return chain.id == "A"

        pdb = CHAIN_A_PDB
        io_writer = PDBIO()
        io_writer.set_structure(s)
        io_writer.save(pdb, ChainASelect())

    print(f"[scan] Loading {pdb}")
    chain, aa_residues, fmn_residue = get_chain_a_residues(pdb, "A")

    # Load the active-site JSON to get the second-shell candidate list
    results_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"
    )
    json_path = os.path.join(results_dir, "nema_active_site.json")
    with open(json_path) as f:
        active_site = json.load(f)
    candidates_meta = active_site["second_shell_candidates"]
    candidate_ids = {c["residue_id"] for c in candidates_meta}
    print(f"[scan] {len(candidate_ids)} second-shell candidates from active-site JSON")

    # Map back to biopython residue objects
    candidates = [r for r in aa_residues if r.id[1] in candidate_ids]
    print(f"[scan] Found biopython residue objects for {len(candidates)} of them")

    # Compute contact-disruption scores
    print(f"[scan] Computing contact-disruption scores (cutoff={CONTACT_CUTOFF_A} A) ...")
    results = scan_candidates(pdb, candidates, aa_residues)
    print(f"[scan] Computed scores for {len(results)} candidates")
    for r in results:
        print(f"[scan]   {r['residue_name']:3s} {r['residue_id']:4d}  "
              f"score={r['contact_disruption_score']:.3f}  band={r['score_band']}")

    # Write CSV
    csv_path = os.path.join(results_dir, "nema_contact_disruption.csv")
    write_csv(results, csv_path)
    print(f"[scan] Wrote {csv_path}")

    # Write figure
    fig_path = os.path.join(results_dir, "nema_contact_disruption.png")
    plot_scores(results, fig_path)
    print(f"[scan] Wrote {fig_path}")

    # Top candidate
    if results:
        top = results[0]
        print(f"[scan] Top candidate (most tolerant): "
              f"{top['residue_name']}{top['residue_id']} "
              f"(score={top['contact_disruption_score']:.3f}, band={top['score_band']})")
        if top["score_band"] == "costly":
            print(f"[scan] WARNING: top candidate is in the 'costly' band "
                  f"(score > {SCORE_COSTLY}). This is a NEGATIVE result for "
                  f"the engineering hypothesis. See manuscript §5.")


if __name__ == "__main__":
    main()
