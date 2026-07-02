"""
prepare_docking.py — Prepares the WT and ILE328ALA mutant NemA receptor
.pdbqt files and the CrO4(2-) ligand .pdbqt for AutoDock Vina docking.

Outputs (in simulations/structural_pipeline/data/):
  - 8BPQ_chainA.pdb              (input, already exists)
  - 8BPQ_chainA_ILE328ALA.pdb    (in-silico mutant, Ile328 -> Ala)
  - 8BPQ_chainA.pdbqt            (receptor WT, for Vina)
  - 8BPQ_chainA_ILE328ALA.pdbqt  (receptor mutant, for Vina)
  - chromate.pdb                 (CrO4(2-) 3D, from rdkit)
  - chromate.pdbqt               (ligand, for Vina)

Vina treats rigid-receptor .pdbqt files as:
  ATOM      1  N   ALA A   1      27.360  24.470  10.000  1.00  0.00    +0.000 N
  Columns:  1-6 record, 7-11 serial, 13-16 atom name, 17 alt, 18-20 resname,
            22 chain, 23-26 resid, 31-38 x, 39-46 y, 47-54 z, 55-60 occ,
            61-66 temp, 77-78 charge, 79-80 atom_type

Atom type is the standard Vina AD4 type inferred from the element:
  C -> C, N -> N, O -> O, S -> S, P -> P, H -> HD,
  most others -> A (acceptor). This is a coarse but functional
  typing for rigid-receptor docking. For higher fidelity, use
  MGLTools' prepare_receptor4.py or ADFRsuite, neither of which
  is available on this machine.

The ILE328 -> ALA in-silico mutation is performed by:
  1. Trimming ILE328's side chain to CB (the Ala-like form)
  2. Renaming the residue to ALA
  This is a real biopython residue replacement, not a SMILES
  rebuild, so the backbone and CB geometry are preserved exactly.
"""

import os
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import iter_heavy_atoms, side_chain_trim_to_ala
from analyze_active_site import CHAIN_A_PDB, FMN_POCKET_CUTOFF
from fetch_nemA_structure import PDB_ID, fetch_8bpq

# Element -> Vina AD4 atom type (coarse but functional for rigid-receptor)
ELEMENT_TO_VINA_TYPE = {
    "C": "C",
    "N": "N",
    "O": "O",
    "S": "S",
    "P": "P",
    "H": "HD",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    "FE": "Fe",
    "ZN": "Zn",
    "MG": "Mg",
    "MN": "Mn",
    "CU": "Cu",
    "CR": "Cr",  # for the chromate Cr atom in the ligand
}


def pdb_to_pdbqt_receptor(pdb_path, pdbqt_path, strip_waters=True):
    """Convert a receptor .pdb to a Vina-compatible .pdbqt.

    For each ATOM/HETATM record (excluding waters), write the
    standard PDB columns 1-66 plus the Vina atom type in columns
    71-76 (right-justified, padded to 6 chars). The standard pdbqt
    format reserves columns 67-70 for the partial charge, which we
    set to 0.000 since this script does not compute Gasteiger
    charges.

    Water HOH residues and the END/TER records are stripped to
    produce a clean protein+cofactor-only receptor that Vina can
    parse without confusion.
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    out_lines = []
    n_written = 0
    for line in lines:
        record = line[:6].strip()
        # Skip water HOH and non-cofactor HETATM
        if strip_waters and record == "HETATM":
            resname = line[17:20].strip().upper()
            if resname == "HOH":
                continue
        if record not in ("ATOM", "HETATM"):
            # Drop TER, END, CONECT, MASTER, etc.
            continue
        # Parse element from columns 77-78 (if present), else from atom name
        element = ""
        if len(line) >= 78:
            element = line[76:78].strip().upper()
        if not element:
            atom_name = line[12:16].strip()
            element = atom_name[0] if atom_name else "C"
        vina_type = ELEMENT_TO_VINA_TYPE.get(element, "A")

        # Build the pdbqt line in the canonical Vina format.
        # Standard pdbqt: 1-6 record, 7-11 serial, 12 space,
        # 13-16 name, 17 altloc, 18-20 resname, 21 space, 22 chain,
        # 23-26 resid, 27 icode, 28-30 spaces, 31-38 x, 39-46 y,
        # 47-54 z, 55-60 occ, 61-66 temp, 67-70 spaces,
        # 71-76 atom type (right-justified, 6 chars).
        try:
            serial = int(line[6:11])
            name = line[12:16]
            altloc = line[16] if len(line) > 16 else " "
            resname = line[17:20]
            chain = line[21] if len(line) > 21 else "A"
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occ = float(line[54:60]) if len(line) >= 60 and line[54:60].strip() else 1.0
            tempf = float(line[60:66]) if len(line) >= 66 and line[60:66].strip() else 0.0
        except ValueError:
            continue

        new_line = (
            f"{record:<6s}"        # 1-6
            f"{serial:>5d} "        # 7-11 serial + space
            f"{name:<4s}"           # 12-16 atom name
            f"{altloc:<1s}"         # 17 altloc
            f"{resname:<3s} "       # 18-20 resname + space
            f"{chain:<1s}"          # 22 chain
            f"{resid:>4d}"          # 23-26 resid
            f"    "                 # 27-30
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"   # 31-54 xyz
            f"{occ:>6.2f}{tempf:>6.2f}"      # 55-66
            f"    "                 # 67-70
            f"{+0.000:>6.3f}"       # 71-76 partial charge (set to 0)
            f" {vina_type:<2s}"     # 77-79 atom type
            f"\n"                   # 80 newline
        )
        out_lines.append(new_line)
        n_written += 1

    with open(pdbqt_path, "w") as f:
        for line in out_lines:
            f.write(line)
    return n_written


def in_silico_mutate_to_ala(input_pdb, output_pdb, chain_id, resid):
    """Mutate a specific residue to ALA by trimming the side chain
    to CB and renaming the residue. Backbone and CB geometry preserved.
    """
    from Bio.PDB import PDBParser, PDBIO, Select
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(PDB_ID, input_pdb)
    chain = s[0][chain_id]
    target = None
    for r in chain.get_residues():
        if r.id[1] == resid:
            target = r
            break
    if target is None:
        raise ValueError(f"Residue {chain_id}:{resid} not found in {input_pdb}")
    if target.get_resname().strip().upper() == "ALA":
        raise ValueError(f"Residue {chain_id}:{resid} is already ALA")
    if target.get_resname().strip().upper() == "GLY":
        raise ValueError(f"Residue {chain_id}:{resid} is GLY; no side chain to trim")

    keep = side_chain_trim_to_ala(target)
    # Remove atoms not in keep
    atoms_to_remove = [a for a in target.get_atoms() if a.name.strip() not in keep]
    for a in atoms_to_remove:
        target.detach_child(a.name)
    # Rename residue
    target.resname = "ALA"
    print(f"[mutate] {chain_id}:{resid} {input_pdb} -> ALA, removed {len(atoms_to_remove)} side-chain atoms")

    io_writer = PDBIO()
    io_writer.set_structure(s)

    class WholeStructure(Select):
        pass

    io_writer.save(output_pdb, WholeStructure())
    return len(atoms_to_remove)


def build_chromate_pdb(pdb_path, seed=42):
    """Build the CrO4(2-) 3D structure from SMILES and save as .pdb.

    Uses UFF (not MMFF) for inorganic atom types, since MMFF does
    not have parameters for Cr.
    """
    chromate_smiles = "[O-][Cr](=O)(=O)[O-]"
    mol = Chem.MolFromSmiles(chromate_smiles)
    if mol is None:
        raise RuntimeError("rdkit could not parse chromate SMILES")
    mol = Chem.AddHs(mol)

    # Embed with ETKDG; use UFF for inorganic-friendly optimization
    # (MMFF throws on Cr atoms; UFF handles them)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise RuntimeError("rdkit could not embed chromate")
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception as e:
        print(f"[chromate] UFF optimization warning: {e} (geometry may be approximate)")

    # Write as PDB
    pdb_block = Chem.MolToPDBBlock(mol)
    with open(pdb_path, "w") as f:
        f.write(pdb_block)
    print(f"[chromate] Wrote {pdb_path} ({mol.GetNumAtoms()} atoms, 1 conformer, seed={seed})")
    return mol


def main():
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    cif, pdb = fetch_8bpq(data_dir)
    if pdb is None or not os.path.exists(pdb):
        raise RuntimeError("Could not get chain A PDB")

    # 1. WT receptor
    print("[prep] Writing WT receptor .pdbqt ...")
    wt_pdbqt = os.path.join(data_dir, f"{PDB_ID}_chainA.pdbqt")
    n = pdb_to_pdbqt_receptor(pdb, wt_pdbqt)
    print(f"[prep] Wrote {wt_pdbqt} ({n} ATOM/HETATM records)")

    # 2. ILE328ALA mutant
    print("[prep] Building ILE328 -> ALA in-silico mutant ...")
    mut_pdb = os.path.join(data_dir, f"{PDB_ID}_chainA_ILE328ALA.pdb")
    in_silico_mutate_to_ala(pdb, mut_pdb, "A", 328)

    print("[prep] Writing mutant receptor .pdbqt ...")
    mut_pdbqt = os.path.join(data_dir, f"{PDB_ID}_chainA_ILE328ALA.pdbqt")
    n = pdb_to_pdbqt_receptor(mut_pdb, mut_pdbqt)
    print(f"[prep] Wrote {mut_pdbqt} ({n} ATOM/HETATM records)")

    # 3. Chromate ligand
    print("[prep] Building CrO4(2-) 3D structure from SMILES ...")
    chromate_pdb = os.path.join(data_dir, "chromate.pdb")
    build_chromate_pdb(chromate_pdb, seed=42)

    print("[prep] Writing ligand .pdbqt ...")
    chromate_pdbqt = os.path.join(data_dir, "chromate.pdbqt")
    n = pdb_to_pdbqt_receptor(chromate_pdb, chromate_pdbqt, strip_waters=False)
    # Add TORSDOF marker for the ligand (rigid, 0 rotatable bonds).
    # Standard pdbqt practice: TORSDOF appears after ENDROOT.
    # For a rigid ligand we wrap everything in ROOT/ENDROOT and
    # append TORSDOF 0.
    with open(chromate_pdbqt) as f:
        body = f.read()
    body = "ROOT\n" + body + "ENDROOT\nTORSDOF 0\n"
    with open(chromate_pdbqt, "w") as f:
        f.write(body)
    print(f"[prep] Wrote {chromate_pdbqt} ({n} ATOM/HETATM records)")


if __name__ == "__main__":
    main()
