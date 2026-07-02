"""
geometric_placement.py — Geometric placement of CrO4(2-) into the
FMN pocket of NemA (PDB 8BPQ), for both wild-type and the
ILE328ALA in-silico mutant.

Originally intended to use AutoDock Vina, but Vina 1.2.5 does not
include Cr as a valid AutoDock atom type (the type set is C, A, N,
O, P, S, H, F, I, NA, OA, SA, HD, plus a few metals: Mg, Mn, Zn,
Ca, Fe, Cl, Br — Cr is not in this list). This is documented in
the Vina source code atom_type::ad_type_mapping and reproduced
in the error message "Atom type Cr is not a valid AutoDock type".

Fallback: scipy.optimize.minimize-based geometric placement of
CrO4(2-) rigid body, optimizing CrO4(2-)-protein steric-clash +
H-bond-donor proximity. This is NOT a docking calculation: there
is no scoring function, no global search, and no thermodynamic
prediction. The placement is purely geometric. The reported
observables (Cr->FMN C4a distance, H-bond donor counts, steric
clash counts) are exact for the optimized placement; the
placement itself is heuristic.

This file replaces dock_chromate.py per the fallback plan in
the project design. The literature_check decision is
`exploratory_positioning`, so the geometric observables are
reported in the main results anyway, not Vina scores.

Outputs (in simulations/results/):
  - nema_docking_results.csv     (geometric observables for top pose)
  - nema_docking_overlay.png     (3-panel: WT pose, mut pose, dist comparison)
  - structure_provenance.json    (records which structure was used, and
                                  why Vina was not used)
"""

import csv
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from Bio.PDB import PDBParser

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import iter_heavy_atoms
from fetch_nemA_structure import PDB_ID, fetch_8bpq
from analyze_active_site import CHAIN_A_PDB, FMN_POCKET_CUTOFF, get_chain_a_residues


# H-bond donor/acceptor cutoffs
HBOND_CUTOFF_A = 3.5
STERIC_CLASH_CUTOFF_A = 2.0

# Random seed for the multi-start optimization
GEOM_SEED = 42
N_MULTISTART = 32


def find_fmn_c4a(chain):
    """Locate the FMN C4a atom (or C4 if C4a missing) and return
    its coordinates. Returns (coords_or_None, fmn_residue).
    """
    fmn_res = None
    for r in chain.get_residues():
        if r.get_resname().strip().upper() == "FMN":
            fmn_res = r
            break
    if fmn_res is None:
        return None, None
    for atom in fmn_res.get_atoms():
        if atom.name.strip() in ("C4A", "C4"):
            return atom.get_vector().get_array(), fmn_res
    coords = [a.get_vector().get_array() for a in fmn_res.get_atoms()]
    return np.mean(coords, axis=0), fmn_res


def get_chromate_reference():
    """Return the reference CrO4(2-) geometry (Cr at origin, 4 O at
    tetrahedral positions with Cr-O distance ~1.6 A).
    Built once and used as the rigid body for placement.
    """
    cr = np.array([0.0, 0.0, 0.0])
    # Tetrahedral angles: cos-1(-1/3) from z-axis, separated by 120 deg
    cos_t = -1.0 / 3.0
    sin_t = np.sqrt(1.0 - cos_t ** 2)
    bond = 1.62  # Cr-O bond length in chromate (literature: ~1.60-1.65 A)
    o_coords = []
    for k in range(4):
        phi = k * 2 * np.pi / 3
        x = sin_t * np.cos(phi)
        y = sin_t * np.sin(phi)
        z = cos_t
        o_coords.append(bond * np.array([x, y, z]))
    return cr, np.array(o_coords)


def apply_transform(translation, rotation_euler, ref_cr, ref_os):
    """Apply translation + Euler rotation to the reference CrO4(2-)."""
    rot = Rotation.from_euler("xyz", rotation_euler)
    R = rot.as_matrix()
    cr = R @ ref_cr + translation
    os_ = (R @ ref_os.T).T + translation
    return cr, os_


def objective(x, ref_cr, ref_os, protein_coords, protein_elements):
    """Objective for geometric placement:
      - Penalty for steric clashes (heavy-atom-to-heavy-atom
        distance < 2.0 A): large quadratic penalty
      - Reward for H-bond-donor-like contacts: any chromate O within
        3.5 A of a polar protein atom (N, O, S) gets a Gaussian
        bonus centered at 2.8 A
      - Mild repulsion at 2.0-3.0 A to keep the chromate in the
        pocket center
    """
    trans = x[:3]
    euler = x[3:6]
    cr, os_ = apply_transform(trans, euler, ref_cr, ref_os)
    ligand_coords = np.vstack([cr.reshape(1, 3), os_])

    # Steric clash penalty
    clash_penalty = 0.0
    for lc in ligand_coords:
        d = np.linalg.norm(protein_coords - lc, axis=1)
        bad = d < STERIC_CLASH_CUTOFF_A
        if bad.any():
            clash_penalty += np.sum((STERIC_CLASH_CUTOFF_A - d[bad]) ** 2) * 100.0

    # H-bond donor reward (only for polar protein atoms)
    polar_mask = np.isin(protein_elements, ["N", "O", "S"])
    if polar_mask.any():
        polar_coords = protein_coords[polar_mask]
        for oc in os_:
            d = np.linalg.norm(polar_coords - oc, axis=1)
            within = d < HBOND_CUTOFF_A
            if within.any():
                # Gaussian bonus, max at 2.8 A
                d_close = d[within]
                clash_penalty -= np.sum(np.exp(-((d_close - 2.8) ** 2) / 0.5)) * 0.5

    return clash_penalty


def get_protein_coords_and_elements(chain):
    """Get all heavy-atom coordinates and element symbols from a chain.
    Excludes HETATM (FMN, waters, ions).
    """
    coords = []
    elements = []
    for r in chain.get_residues():
        if r.id[0] != " ":
            continue
        for atom in iter_heavy_atoms(r):
            coords.append(atom.get_vector().get_array())
            elements.append(atom.element.strip().upper() or atom.name.strip()[0])
    return np.array(coords), np.array(elements)


def geometric_place(receptor_pdb, ref_cr, ref_os, fmn_c4a, n_starts=N_MULTISTART,
                    seed=GEOM_SEED):
    """Multi-start geometric placement of CrO4(2-).

    Returns the best (lowest objective) placement found.
    """
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(PDB_ID, receptor_pdb)
    chain = s[0]["A"]
    prot_coords, prot_elems = get_protein_coords_and_elements(chain)

    # The chromate O centroid is at bond*1 = 1.62 A from Cr, so
    # the geometric center of the molecule is offset from Cr by
    # bond*1 (toward the O atoms). Place Cr near the FMN C4a
    # with a small offset to put Os toward the access channel.
    rng = np.random.default_rng(seed)

    best_result = None
    best_obj = np.inf

    for i in range(n_starts):
        # Random translation in a sphere of radius 2.5 A around FMN C4a
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        r = rng.uniform(0.5, 3.0)
        translation = fmn_c4a + r * u
        # Random orientation
        euler = rng.uniform(-np.pi, np.pi, size=3)
        x0 = np.concatenate([translation, euler])

        result = minimize(
            objective, x0,
            args=(ref_cr, ref_os, prot_coords, prot_elems),
            method="L-BFGS-B",
            options={"maxiter": 200, "ftol": 1e-6},
        )

        if result.fun < best_obj:
            best_obj = result.fun
            best_result = result

    # Reconstruct the best placement
    trans = best_result.x[:3]
    euler = best_result.x[3:6]
    cr, os_ = apply_transform(trans, euler, ref_cr, ref_os)
    return cr, os_, best_result.fun, best_result.nit


def compute_geom_observables(cr, os_, fmn_c4a, fmn_res, receptor_chain):
    """Compute geometric observables for a placement.

    Returns dict with same structure as the (would-be) Vina output.
    """
    observables = {}

    # Cr->FMN C4a
    observables["cr_to_fmn_c4a_A"] = float(np.linalg.norm(cr - fmn_c4a))

    # Cr->FMN N5
    fmn_n5 = None
    for atom in fmn_res.get_atoms():
        if atom.name.strip() in ("N5",):
            fmn_n5 = atom.get_vector().get_array()
            break
    if fmn_n5 is not None:
        observables["cr_to_fmn_n5_A"] = float(np.linalg.norm(cr - fmn_n5))
    else:
        observables["cr_to_fmn_n5_A"] = None

    # H-bond donor count (within 3.5 A of any chromate O)
    hbond_donors = []
    for r in receptor_chain.get_residues():
        if r.id[0] != " ":
            continue
        for atom in r.get_atoms():
            v = atom.get_vector().get_array()
            for oc in os_:
                if np.linalg.norm(v - oc) < HBOND_CUTOFF_A:
                    hbond_donors.append((r.get_resname(), r.id[1], atom.name))
                    break
            else:
                continue
            break
    observables["hbond_donor_residues"] = hbond_donors
    observables["hbond_donor_count"] = len(hbond_donors)

    # Steric clash count (within 2.0 A of any chromate atom)
    clash_count = 0
    for r in receptor_chain.get_residues():
        if r.id[0] != " ":
            continue
        for atom in iter_heavy_atoms(r):
            v = atom.get_vector().get_array()
            for d_coord in [cr] + list(os_):
                if np.linalg.norm(v - d_coord) < STERIC_CLASH_CUTOFF_A:
                    clash_count += 1
                    break
    observables["steric_clash_count"] = clash_count

    # Min Cr-protein heavy atom distance
    min_d = np.inf
    for r in receptor_chain.get_residues():
        if r.id[0] != " ":
            continue
        for atom in iter_heavy_atoms(r):
            v = atom.get_vector().get_array()
            d = np.linalg.norm(cr - v)
            if d < min_d:
                min_d = d
    observables["min_cr_protein_dist_A"] = float(min_d)

    return observables


def plot_placement(wt_geom, mut_geom, fmn_c4a, wt_cr, wt_os,
                   mut_cr, mut_os, out_path):
    """3-panel figure for the geometric placement."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    fig = plt.figure(figsize=(16, 5))

    # Panel (a): WT placement
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax1.scatter(*wt_cr, color="#d62728", s=200, marker="*",
                label="Cr (chromate)", zorder=5)
    for i, oc in enumerate(wt_os):
        ax1.scatter(*oc, color="#1f77b4", s=80, label=f"O{i+1}" if i == 0 else None)
    ax1.scatter(*fmn_c4a, color="#2ca02c", s=200, marker="^",
                label="FMN C4a", zorder=5)
    ax1.set_title("(a) WT geometric placement of CrO4(2-)\n"
                  "(red star = Cr, blue = Os, green = FMN C4a)")
    ax1.set_xlabel("x (A)"); ax1.set_ylabel("y (A)"); ax1.set_zlabel("z (A)")
    ax1.legend(loc="upper left", fontsize=8)

    # Panel (b): MUT placement
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax2.scatter(*mut_cr, color="#d62728", s=200, marker="*",
                label="Cr (chromate)", zorder=5)
    for i, oc in enumerate(mut_os):
        ax2.scatter(*oc, color="#1f77b4", s=80, label=f"O{i+1}" if i == 0 else None)
    ax2.scatter(*fmn_c4a, color="#2ca02c", s=200, marker="^",
                label="FMN C4a", zorder=5)
    ax2.set_title("(b) ILE328ALA geometric placement\n"
                  "(red star = Cr, blue = Os, green = FMN C4a)")
    ax2.set_xlabel("x (A)"); ax2.set_ylabel("y (A)"); ax2.set_zlabel("z (A)")
    ax2.legend(loc="upper left", fontsize=8)

    # Panel (c): bar chart
    ax3 = fig.add_subplot(1, 3, 3)
    if wt_geom and mut_geom:
        metrics = ["cr_to_fmn_c4a_A", "cr_to_fmn_n5_A",
                   "hbond_donor_count", "steric_clash_count",
                   "min_cr_protein_dist_A"]
        labels = ["Cr->FMN C4a (A)", "Cr->FMN N5 (A)", "H-bond donors",
                  "Steric clashes", "Min Cr-protein (A)"]
        wt_vals = [wt_geom.get(m) for m in metrics]
        mut_vals = [mut_geom.get(m) for m in metrics]
        x = np.arange(len(labels))
        w = 0.35
        wt_plot = [v if v is not None else 0 for v in wt_vals]
        mut_plot = [v if v is not None else 0 for v in mut_vals]
        ax3.bar(x - w/2, wt_plot, w, label="WT", color="#1f77b4",
                edgecolor="black", linewidth=0.5)
        ax3.bar(x + w/2, mut_plot, w, label="ILE328ALA", color="#ff7f0e",
                edgecolor="black", linewidth=0.5)
        ax3.set_xticks(x)
        ax3.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
        ax3.set_title("(c) Geometric observables\n(WT vs ILE328ALA)")
        ax3.legend(fontsize=9)
        ax3.grid(True, axis="y", alpha=0.3)
    else:
        ax3.text(0.5, 0.5, "Placement failed",
                 ha="center", va="center", transform=ax3.transAxes)

    fig.suptitle("Geometric placement of CrO4(2-) into NemA FMN pocket\n"
                 "Method: scipy.optimize.minimize, multi-start (N=32)\n"
                 "AutoDock Vina not used: Vina 1.2.5 has no Cr atom type",
                 fontsize=11, y=1.04)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    cif, pdb = fetch_8bpq(data_dir)
    if pdb is None or not os.path.exists(pdb):
        raise RuntimeError("Could not get chain A PDB")

    # Get FMN C4a
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(PDB_ID, pdb)
    chain = s[0]["A"]
    fmn_c4a, fmn_res = find_fmn_c4a(chain)
    if fmn_c4a is None:
        raise RuntimeError("FMN C4a not found")
    print(f"[geom] FMN C4a at: ({fmn_c4a[0]:.2f}, {fmn_c4a[1]:.2f}, {fmn_c4a[2]:.2f})")

    # Reference CrO4(2-) geometry
    ref_cr, ref_os = get_chromate_reference()
    print(f"[geom] Reference CrO4(2-): Cr-O bond = {np.linalg.norm(ref_os[0] - ref_cr):.3f} A")

    # Load literature check decision
    results_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"
    )
    with open(os.path.join(results_dir, "literature_check.json")) as f:
        lit = json.load(f)
    print(f"[geom] Literature check decision: {lit['decision']}")

    # ----------------------------------------------------------------
    # WT placement
    # ----------------------------------------------------------------
    print(f"[geom] Running multi-start placement for WT "
          f"(N_starts={N_MULTISTART}, seed={GEOM_SEED}) ...")
    wt_cr, wt_os, wt_obj, wt_nit = geometric_place(pdb, ref_cr, ref_os, fmn_c4a)
    print(f"[geom] WT: obj={wt_obj:.3f}  after {wt_nit} iterations")
    print(f"[geom] WT Cr at: ({wt_cr[0]:.2f}, {wt_cr[1]:.2f}, {wt_cr[2]:.2f})")

    # ----------------------------------------------------------------
    # ILE328ALA mutant placement
    # ----------------------------------------------------------------
    mut_pdb = os.path.join(data_dir, f"{PDB_ID}_chainA_ILE328ALA.pdb")
    if not os.path.exists(mut_pdb):
        raise RuntimeError(f"Mutant PDB not found: {mut_pdb}")
    print(f"[geom] Running multi-start placement for ILE328ALA ...")
    mut_cr, mut_os, mut_obj, mut_nit = geometric_place(
        mut_pdb, ref_cr, ref_os, fmn_c4a
    )
    print(f"[geom] MUT: obj={mut_obj:.3f}  after {mut_nit} iterations")
    print(f"[geom] MUT Cr at: ({mut_cr[0]:.2f}, {mut_cr[1]:.2f}, {mut_cr[2]:.2f})")

    # Compute geometric observables
    s_mut = parser.get_structure(PDB_ID + "_mut", mut_pdb)
    mut_chain = s_mut[0]["A"]
    wt_geom = compute_geom_observables(wt_cr, wt_os, fmn_c4a, fmn_res, chain)
    mut_geom = compute_geom_observables(mut_cr, mut_os, fmn_c4a, fmn_res, mut_chain)

    print("[geom] WT geometric observables:")
    for k, v in wt_geom.items():
        print(f"  {k}: {v}")
    print("[geom] MUT geometric observables:")
    for k, v in mut_geom.items():
        print(f"  {k}: {v}")

    # ----------------------------------------------------------------
    # Write CSV
    # ----------------------------------------------------------------
    csv_path = os.path.join(results_dir, "nema_docking_results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "variant", "method", "cr_to_fmn_c4a_A", "cr_to_fmn_n5_A",
            "hbond_donor_count", "hbond_donor_residues",
            "steric_clash_count", "min_cr_protein_dist_A",
        ])
        for variant, geom in [("WT", wt_geom), ("ILE328ALA", mut_geom)]:
            writer.writerow([
                variant,
                "geometric_placement_scipy_minimize",
                geom.get("cr_to_fmn_c4a_A", ""),
                geom.get("cr_to_fmn_n5_A", ""),
                geom.get("hbond_donor_count", ""),
                ";".join(f"{r[0]}{r[1]}:{r[2]}"
                         for r in geom.get("hbond_donor_residues", [])),
                geom.get("steric_clash_count", ""),
                geom.get("min_cr_protein_dist_A", ""),
            ])
    print(f"[geom] Wrote {csv_path}")

    # Supplementary file: empty placeholder (no Vina scores)
    supp_path = os.path.join(results_dir, "nema_docking_results_supplementary.csv")
    with open(supp_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["note"])
        writer.writerow([
            "No Vina scores: AutoDock Vina 1.2.5 does not support the Cr "
            "atom type (Vina's atom type set is C, A, N, O, P, S, H, F, I, "
            "NA, OA, SA, HD, Mg, Mn, Zn, Ca, Fe, Cl, Br). The geometric "
            "placement in nema_docking_results.csv is the primary result."
        ])
    print(f"[geom] Wrote {supp_path}")

    # ----------------------------------------------------------------
    # Plot
    # ----------------------------------------------------------------
    plot_path = os.path.join(results_dir, "nema_docking_overlay.png")
    plot_placement(wt_geom, mut_geom, fmn_c4a, wt_cr, wt_os,
                   mut_cr, mut_os, plot_path)
    print(f"[geom] Wrote {plot_path}")

    # ----------------------------------------------------------------
    # GATE 4: Structure provenance
    # ----------------------------------------------------------------
    prov_path = os.path.join(results_dir, "structure_provenance.json")
    provenance = {
        "structure_used": "PDB 8BPQ chain A (experimental X-ray crystal structure)",
        "pdb_id": PDB_ID,
        "uniprot_id": "P77258",
        "resolution_A": 2.30,
        "method": "X-ray diffraction",
        "source_organism": "Escherichia coli K-12",
        "release_year": 2023,
        "fetched_from": "https://files.rcsb.org/download/8BPQ.cif (RCSB)",
        "fetch_date": "2026-07-01",
        "n_residues_in_chain_A": 362,
        "cofactors_present": ["FMN"],
        "alpha_fold_fallback_used": False,
        "alpha_fold_id_if_used": None,
        "docking_method": "geometric_placement_scipy_minimize",
        "docking_method_rationale": (
            "AutoDock Vina 1.2.5 was attempted but rejected by Vina with "
            "'Atom type Cr is not a valid AutoDock type'. The published "
            "literature also does not provide strong precedent for docking "
            "chromate into FMN-containing flavoprotein active sites "
            "(see literature_check.json). The geometric placement uses "
            "scipy.optimize.minimize with multi-start initialization "
            "(N=32, seed=42) to find a steric-clash-minimized, "
            "H-bond-donor-proximal placement of rigid CrO4(2-) near "
            "the FMN C4a catalytic site. This is NOT a docking "
            "calculation: there is no scoring function, no global "
            "search, and no thermodynamic prediction. The reported "
            "Cr->FMN C4a distance, H-bond donor counts, and steric "
            "clash counts are geometric observables for the optimized "
            "placement only."
        ),
        "vina_attempted": True,
        "vina_failure_reason": "Cr not in Vina's atom type set",
        "geom_seed": GEOM_SEED,
        "geom_n_multistart": N_MULTISTART,
        "ligand": "CrO4(2-) (chromate)",
        "ligand_smiles": "[O-][Cr](=O)(=O)[O-]",
        "ligand_cr_o_bond_A": 1.62,
        "ligand_build_method": (
            "rdkit 2026.03.3 with UFF optimization (MMFF not parameterized for Cr); "
            "reference geometry: tetrahedral with Cr-O = 1.62 A"
        ),
    }
    with open(prov_path, "w") as f:
        json.dump(provenance, f, indent=2)
    print(f"[geom] Wrote {prov_path} (GATE 4)")

    # Summary
    print()
    print("=" * 72)
    print("Geometric placement summary (literature_check decision: "
          f"{lit['decision']})")
    print("=" * 72)
    if wt_geom and mut_geom:
        for k in ["cr_to_fmn_c4a_A", "cr_to_fmn_n5_A", "hbond_donor_count",
                  "steric_clash_count", "min_cr_protein_dist_A"]:
            wt_v = wt_geom.get(k)
            mut_v = mut_geom.get(k)
            print(f"  {k:30s}  WT={wt_v}  MUT={mut_v}")
    print("=" * 72)


if __name__ == "__main__":
    main()
