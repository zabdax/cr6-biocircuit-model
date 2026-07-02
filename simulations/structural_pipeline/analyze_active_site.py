"""
analyze_active_site.py — Identifies the FMN binding pocket of E. coli
NemA (PDB 8BPQ, chain A) and classifies pocket residues by their
relationship to the FMN cofactor and the UniProt-annotated catalytic
machinery.

Outputs:
  - simulations/results/nema_active_site.json   (pocket composition)
  - simulations/results/nema_active_site.png   (3-panel figure: full
    chain A trace, zoomed FMN pocket with catalytic residues, second-
    shell residue list with side-chain size classification)

Definitions used in this analysis:
  FMN binding pocket:    all residues with any heavy atom within 5.0 A
                        of any FMN heavy atom
  FMN-binding:           residues with any heavy atom within 3.5 A
                        of any FMN heavy atom (direct contacts)
  First-shell catalytic: residues annotated by UniProt as
                        FMN-binding (101, 234, 324, 325) or as
                        the proton donor (187), per P77258
  Second-shell:          pocket residues (5.0 A) with large
                        aromatic / hydrophobic side chains
                        (TRP, PHE, TYR, ILE, LEU, MET) that are
                        not first-shell catalytic. These are the
                        candidate set for the NemA*2+ engineering
                        hypothesis (Mowafy 2010 OYE precedent).

References:
  - UniProt P77258: https://rest.uniprot.org/uniprotkb/P77258
  - PDB 8BPQ:        https://www.rcsb.org/structure/8BPQ
  - Mowafy 2010:     Protein Sci 19(7):1283-1296 (OYE active-site
                     mutagenesis precedent for second-shell pocket
                     engineering)
"""

import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import iter_heavy_atoms, iter_residue_heavy_coords

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fetch_nemA_structure import PDB_ID, fetch_8bpq


# UniProt P77258 catalytic-residue annotations
UNIPROT_FMN_BINDING_RESIDUES = {101, 234, 324, 325}   # FMN contact
UNIPROT_PROTON_DONOR_RESIDUES = {187}                  # proton donor
UNIPROT_CATALYTIC = UNIPROT_FMN_BINDING_RESIDUES | UNIPROT_PROTON_DONOR_RESIDUES

# Side-chain categories
AROMATIC_HYDROPHOBIC = {"TRP", "PHE", "TYR", "ILE", "LEU", "MET", "VAL"}
CHARGED_POSITIVE = {"ARG", "LYS", "HIS"}
CHARGED_NEGATIVE = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "CYS"}

# Contact cutoffs (A)
FMN_POCKET_CUTOFF = 5.0
FMN_BINDING_CUTOFF = 3.5

# Structure file locations
CHAIN_A_PDB = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", f"{PDB_ID}_chainA.pdb"
)


def get_chain_a_residues(structure_path, chain_id="A"):
    """Load a PDB file via biopython and return the residue iterator
    for the specified chain, excluding waters and other HETATM.

    Returns (chain, all_residues, fmn_residue_or_None).
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(PDB_ID, structure_path)
    model = s[0]
    chain = model[chain_id]

    aa_residues = [r for r in chain.get_residues() if r.id[0] == " "]
    het_residues = [r for r in chain.get_residues() if r.id[0] != " "]

    fmn_residue = None
    for r in het_residues:
        if r.get_resname().strip().upper() == "FMN":
            fmn_residue = r
            break

    return chain, aa_residues, fmn_residue


def pocket_residues_from_fmn(aa_residues, fmn_residue, cutoff=FMN_POCKET_CUTOFF):
    """Return the list of (residue, min_dist_to_FMN) for AA residues
    with any heavy atom within `cutoff` of any FMN heavy atom.
    """
    fmn_coords = [c for _, c in iter_residue_heavy_coords(fmn_residue)]
    fmn_arr = np.array(fmn_coords)

    out = []
    for r in aa_residues:
        r_coords = [c for _, c in iter_residue_heavy_coords(r)]
        if not r_coords:
            continue
        r_arr = np.array(r_coords)
        diff = r_arr[:, None, :] - fmn_arr[None, :, :]
        d = np.sqrt(np.sum(diff * diff, axis=2))
        dmin = float(d.min())
        if dmin <= cutoff:
            out.append((r, dmin))
    out.sort(key=lambda x: x[1])
    return out


def classify_residue(res):
    """Classify a residue by side-chain category."""
    rn = res.get_resname().strip().upper()
    if rn in AROMATIC_HYDROPHOBIC:
        return "aromatic_hydrophobic"
    if rn in CHARGED_POSITIVE:
        return "charged_positive"
    if rn in CHARGED_NEGATIVE:
        return "charged_negative"
    if rn in POLAR:
        return "polar"
    return "other"


def annotate_pocket(pocket):
    """For each (residue, dist) entry, return a dict with full annotation.

    Fields:
      residue_id, residue_name, chain, dist_to_FMN_A,
      is_fmn_binding, is_first_shell_catalytic, is_second_shell,
      side_chain_class

    Definitions:
      is_fmn_binding:           min heavy-atom distance to FMN <= 3.5 A
      is_first_shell_catalytic: in UniProt-annotated catalytic set
                                (FMN-binding residues 101/234/324/325 or
                                proton donor 187) OR within 3.5 A of FMN
                                (direct FMN contact)
      is_second_shell_candidate: pocket residue (within 5.0 A) that is
                                aromatic/hydrophobic AND is NOT a
                                first-shell catalytic residue. These
                                are the engineering candidates the
                                Mowafy 2010 OYE precedent targets
                                (second-shell pocket enlargement).
    """
    out = []
    for r, d in pocket:
        rn = r.get_resname().strip().upper()
        rid = r.id[1]
        sc = classify_residue(r)
        is_fmn_binding = d <= FMN_BINDING_CUTOFF
        # First-shell = UniProt-annotated OR direct FMN contact (<=3.5 A).
        # This excludes TYR 352 and LEU 26 (which contact FMN at ~3.2-3.4 A
        # but are not UniProt-annotated) from the second-shell set, since
        # they are direct FMN contacts rather than flanking residues.
        is_first_shell = (rid in UNIPROT_CATALYTIC) or is_fmn_binding
        is_second_shell = (
            sc == "aromatic_hydrophobic" and not is_first_shell
        )
        out.append({
            "residue_id": int(rid),
            "residue_name": rn,
            "chain": r.get_parent().id,
            "dist_to_FMN_A": round(d, 2),
            "is_fmn_binding": bool(is_fmn_binding),
            "is_first_shell_catalytic": bool(is_first_shell),
            "is_second_shell_candidate": bool(is_second_shell),
            "side_chain_class": sc,
        })
    return out


def write_json(annotated, out_path):
    """Write annotated pocket to JSON with summary header."""
    fmn_b = [r for r in annotated if r["is_fmn_binding"]]
    first = [r for r in annotated if r["is_first_shell_catalytic"]]
    second = [r for r in annotated if r["is_second_shell_candidate"]]

    payload = {
        "pdb_id": PDB_ID,
        "source_structure": "PDB 8BPQ chain A (2.30 A X-ray, UniProt P77258)",
        "contact_cutoffs_A": {
            "pocket_radius": FMN_POCKET_CUTOFF,
            "fmn_binding_radius": FMN_BINDING_CUTOFF,
        },
        "uniprot_annotations": {
            "FMN_binding_residues": sorted(UNIPROT_FMN_BINDING_RESIDUES),
            "proton_donor_residues": sorted(UNIPROT_PROTON_DONOR_RESIDUES),
            "all_catalytic": sorted(UNIPROT_CATALYTIC),
        },
        "summary": {
            "n_pocket_residues": len(annotated),
            "n_fmn_binding": len(fmn_b),
            "n_first_shell_catalytic": len(first),
            "n_second_shell_candidates": len(second),
        },
        "pocket_residues": annotated,
        "first_shell_catalytic_in_pocket": first,
        "second_shell_candidates": second,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    return payload


def plot_pocket(payload, structure_path, out_path):
    """3-panel figure:
      (a) full chain A C-alpha trace in light grey, FMN-binding
          residues in red, first-shell catalytic in orange, second-
          shell candidates in blue
      (b) zoomed view of the FMN pocket (same color scheme)
      (c) bar chart of second-shell candidates, x = residue_id,
          y = distance to FMN
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(PDB_ID, structure_path)
    model = s[0]
    chain = model["A"]

    # Extract C-alpha coordinates
    ca_coords = []
    ca_residues = []
    for r in chain.get_residues():
        if r.id[0] != " ":
            continue
        if "CA" in r:
            ca_coords.append(r["CA"].get_vector().get_array())
            ca_residues.append(r)
    ca_arr = np.array(ca_coords)

    # FMN centroid
    fmn_res = None
    for r in chain.get_residues():
        if r.get_resname().strip().upper() == "FMN":
            fmn_res = r
            break
    fmn_coords = [a.get_vector().get_array() for a in fmn_res.get_atoms()]
    fmn_centroid = np.mean(fmn_coords, axis=0)
    fmn_C4a = None
    for atom in fmn_res.get_atoms():
        if atom.name.strip() in ("C4A", "C4"):
            fmn_C4a = atom.get_vector().get_array()
            break
    if fmn_C4a is None:
        # Approximate C4a as the centroid if exact atom missing
        fmn_C4a = fmn_centroid

    # Categorize CA residues
    fmn_b_ids = {r["residue_id"] for r in payload["pocket_residues"]
                 if r["is_fmn_binding"]}
    first_ids = {r["residue_id"] for r in payload["first_shell_catalytic_in_pocket"]}
    second_ids = {r["residue_id"] for r in payload["second_shell_candidates"]}
    pocket_ids = {r["residue_id"] for r in payload["pocket_residues"]}

    def get_color(r):
        rid = r.id[1]
        if rid in fmn_b_ids:
            return "#d62728"  # red: direct FMN contact (<=3.5 A)
        if rid in second_ids:
            return "#1f77b4"  # blue: second-shell aromatic/hydrophobic
        if rid in pocket_ids:
            return "#ff7f0e"  # orange: other pocket residue (charged/polar)
        return None  # non-pocket: do not highlight (keep grey trace)

    fig = plt.figure(figsize=(16, 5))
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax3 = fig.add_subplot(1, 3, 3)

    # Panel (a): full chain trace
    ax1.plot(ca_arr[:, 0], ca_arr[:, 1], ca_arr[:, 2], "-",
             color="#999999", lw=0.5, alpha=0.6)
    for i, r in enumerate(ca_residues):
        c = get_color(r)
        if c is not None:
            ax1.scatter(*ca_arr[i], color=c, s=30, zorder=3)
    ax1.scatter(*fmn_centroid, color="#2ca02c", s=120, marker="*",
                label="FMN centroid", zorder=5)
    ax1.set_title("(a) Full chain A — colored: FMN binding (red),\n"
                  "other pocket residues (orange), 2nd-shell candidates (blue)")
    ax1.set_xlabel("x (A)"); ax1.set_ylabel("y (A)"); ax1.set_zlabel("z (A)")
    ax1.legend(loc="upper left", fontsize=8)

    # Panel (b): zoomed view (residues within 15 A of FMN C4a)
    zoom_radius = 15.0
    dist_to_c4a = np.linalg.norm(ca_arr - fmn_C4a, axis=1)
    in_zoom = dist_to_c4a < zoom_radius
    ax2.plot(ca_arr[in_zoom, 0], ca_arr[in_zoom, 1], ca_arr[in_zoom, 2],
             "-", color="#999999", lw=1.0, alpha=0.8)
    for i, r in enumerate(ca_residues):
        if not in_zoom[i]:
            continue
        c = get_color(r)
        if c is not None:
            ax2.scatter(*ca_arr[i], color=c, s=80, zorder=3, edgecolor="black",
                        linewidth=0.5)
    ax2.scatter(*fmn_C4a, color="#2ca02c", s=200, marker="*",
                label="FMN C4a (catalytic)", zorder=5)
    # FMN heavy atoms in green dots
    for a in fmn_res.get_atoms():
        v = a.get_vector().get_array()
        ax2.scatter(*v, color="#2ca02c", s=15, alpha=0.5)
    ax2.set_title(f"(b) FMN pocket zoom (within {zoom_radius:.0f} A of FMN C4a)")
    ax2.set_xlabel("x (A)"); ax2.set_ylabel("y (A)"); ax2.set_zlabel("z (A)")
    ax2.legend(loc="upper left", fontsize=8)

    # Panel (c): bar chart of second-shell candidates
    second = payload["second_shell_candidates"]
    if second:
        ids = [r["residue_id"] for r in second]
        dists = [r["dist_to_FMN_A"] for r in second]
        names = [r["residue_name"] for r in second]
        labels = [f"{n}{i}" for n, i in zip(names, ids)]
        # All second-shell candidates are by definition NOT within the
        # FMN-binding cutoff (3.5 A), so all bars are blue.
        ax3.barh(labels, dists, color="#1f77b4", edgecolor="black", linewidth=0.5)
        ax3.axvline(FMN_BINDING_CUTOFF, color="orange", linestyle="--",
                    lw=1.5, label=f"FMN-binding cutoff ({FMN_BINDING_CUTOFF} A)")
        ax3.set_xlabel("min heavy-atom distance to FMN (A)")
        ax3.set_title("(c) Second-shell aromatic/hydrophobic residues\n"
                      "flanking the FMN pocket (all > 3.5 A from FMN)")
        ax3.invert_yaxis()
        ax3.grid(True, axis="x", alpha=0.3)
        ax3.legend(loc="lower right", fontsize=8)
    else:
        ax3.text(0.5, 0.5, "No second-shell candidates found",
                 ha="center", va="center", transform=ax3.transAxes)
        ax3.set_title("(c) Second-shell candidates")

    fig.suptitle(
        f"NemA (PDB {PDB_ID}) active-site analysis — FMN pocket composition",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    # Ensure structure is fetched
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    cif, pdb = fetch_8bpq(data_dir)
    if pdb is None or not os.path.exists(pdb):
        # Fall back to extracting from cif
        from Bio.PDB import MMCIFParser, PDBIO, Select
        parser = MMCIFParser(QUIET=True)
        s = parser.get_structure(PDB_ID, cif)

        class ChainASelect(Select):
            def accept_chain(self, chain):
                return chain.id == "A"

        tmp_pdb = CHAIN_A_PDB
        io_writer = PDBIO()
        io_writer.set_structure(s)
        io_writer.save(tmp_pdb, ChainASelect())
        pdb = tmp_pdb

    print(f"[analyze] Loading {pdb}")
    chain, aa_residues, fmn_residue = get_chain_a_residues(pdb, "A")
    if fmn_residue is None:
        raise RuntimeError("FMN cofactor not found in chain A of 8BPQ")
    print(f"[analyze] Chain A: {len(aa_residues)} AA residues")
    print(f"[analyze] FMN cofactor: {fmn_residue.get_resname()} "
          f"with {len(list(fmn_residue.get_atoms()))} atoms")

    # Compute pocket
    pocket = pocket_residues_from_fmn(aa_residues, fmn_residue,
                                      cutoff=FMN_POCKET_CUTOFF)
    print(f"[analyze] Pocket: {len(pocket)} residues within "
          f"{FMN_POCKET_CUTOFF} A of FMN")
    annotated = annotate_pocket(pocket)
    n_first = sum(1 for r in annotated if r["is_first_shell_catalytic"])
    n_second = sum(1 for r in annotated if r["is_second_shell_candidate"])
    print(f"[analyze]   FMN-binding (<=3.5 A): "
          f"{sum(1 for r in annotated if r['is_fmn_binding'])}")
    print(f"[analyze]   First-shell catalytic: {n_first}")
    print(f"[analyze]   Second-shell candidates: {n_second}")
    print("[analyze] Second-shell candidates:")
    for r in annotated:
        if r["is_second_shell_candidate"]:
            print(f"[analyze]   {r['residue_name']:3s} {r['residue_id']:4d}  "
                  f"dist={r['dist_to_FMN_A']:.2f} A  "
                  f"class={r['side_chain_class']}")

    # Write JSON
    results_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"
    )
    os.makedirs(results_dir, exist_ok=True)
    json_path = os.path.join(results_dir, "nema_active_site.json")
    payload = write_json(annotated, json_path)
    print(f"[analyze] Wrote {json_path}")

    # Write figure
    fig_path = os.path.join(results_dir, "nema_active_site.png")
    plot_pocket(payload, pdb, fig_path)
    print(f"[analyze] Wrote {fig_path}")


if __name__ == "__main__":
    main()
