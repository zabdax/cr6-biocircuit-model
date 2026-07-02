"""
utils.py — Shared helpers for the SJWP 2026 NemA*2+ structural pipeline.

These helpers are deliberately minimal: they wrap biopython idioms
that recur across multiple pipeline stages and apply the project's
heavy-atom / contact-cutoff conventions consistently.
"""

import numpy as np


# Heavy-atom names (everything except H). Used throughout the
# pipeline to keep distance and contact calculations consistent
# (the standard biopython .get_atoms() iterator includes H if
# present in the PDB, which would bias contact counts).
HEAVY_ATOM_NAMES = {
    "C", "N", "O", "S", "P", "SE", "FE", "ZN", "MG", "MN", "CU", "NI",
    "CO", "MO", "W", "F", "CL", "BR", "I", "CR",  # CR for chromate Cr
}


def iter_heavy_atoms(residue):
    """Yield only heavy (non-hydrogen) atoms from a biopython Residue."""
    for atom in residue.get_atoms():
        if atom.element.strip().upper() in HEAVY_ATOM_NAMES or \
           atom.name.strip().startswith(("C", "N", "O", "S", "P")):
            # Take any atom whose element symbol or first-letter-of-name
            # is a heavy element. This is robust to structures with
            # partially populated element fields.
            yield atom


def iter_residue_heavy_coords(residue):
    """Yield (atom, np.array([x,y,z])) for each heavy atom in a residue."""
    for atom in iter_heavy_atoms(residue):
        yield atom, atom.get_vector().get_array()


def min_distance_between_residues(res_a, res_b):
    """Minimum heavy-atom-to-heavy-atom distance between two residues,
    in Angstroms. Returns np.inf if either residue has no heavy atoms.
    """
    coords_a = [c for _, c in iter_residue_heavy_coords(res_a)]
    coords_b = [c for _, c in iter_residue_heavy_coords(res_b)]
    if not coords_a or not coords_b:
        return np.inf
    a = np.array(coords_a)
    b = np.array(coords_b)
    # Pairwise distance matrix
    diff = a[:, None, :] - b[None, :, :]
    d = np.sqrt(np.sum(diff * diff, axis=2))
    return float(d.min())


def side_chain_trim_to_ala(residue):
    """Return a list of atom names to KEEP if we trim this residue's
    side chain to a methyl (alanine-like). N, CA, C, O, plus CB.
    For Gly, only the backbone remains (no CB).
    """
    resname = residue.get_resname().strip().upper()
    keep = {"N", "CA", "C", "O"}
    if resname != "GLY":
        keep.add("CB")
    return keep
