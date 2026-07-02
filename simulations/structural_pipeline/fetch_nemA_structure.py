"""
fetch_nemA_structure.py — Fetches the E. coli NemA crystal structure
(PDB 8BPQ, 2.30 A X-ray, UniProt P77258) for downstream structural
analysis in the SJWP 2026 NemA*2+ evaluation pipeline.

Output:
  - simulations/structural_pipeline/data/8BPQ.cif       (mmCIF, uncompressed)
  - simulations/structural_pipeline/data/8BPQ_chainA.pdb  (chain A, .pdb)

The mmCIF file is fetched directly from RCSB via urllib (stdlib only),
since the pdb-database skill's helper script requires the `uv` runner
which is not installed on this Windows machine. The fetch URL pattern
follows the same shard-based layout used by the skill's
download_coordinate_files.py script:

  https://files.rcsb.org/download/<PDBID>.cif

8BPQ is a 2.30 A X-ray structure of E. coli NemA. UniProt P77258
cross-references 8BPQ (2.30 A) and 8BPP (3.10 A). 8BPQ is the higher-
resolution structure and is the canonical choice for active-site
analysis. Note: 1OFW (which appears in some prior literature as a
NemA reference) is NOT NemA — it is the nine-heme cytochrome c from
Desulfovibrio desulfuricans. This script deliberately targets 8BPQ.

References:
  - UniProt P77258: https://rest.uniprot.org/uniprotkb/P77258
  - RCSB PDB 8BPQ:  https://www.rcsb.org/structure/8BPQ
  - RCSB usage policy: https://www.rcsb.org/pages/usage-policy
"""

import gzip
import os
import sys
import urllib.request

# The fetch URL pattern follows the RCSB shard layout.
# For 8BPQ, the sharded path under files.rcsb.org/download/ is:
#   8BPQ.cif
PDB_ID = "8BPQ"
RCSB_DOWNLOAD_URL = f"https://files.rcsb.org/download/{PDB_ID}.cif"
# Mirror with .cif.gz in case the uncompressed file is unavailable
RCSB_DOWNLOAD_URL_GZ = f"https://files.rcsb.org/download/{PDB_ID}.cif.gz"

# Fallback to the older pdbj mirror if files.rcsb.org is unreachable
PDBJ_FALLBACK_URL = f"https://pdbj.org/rest/download/{PDB_ID}.cif"


def fetch_url(url, timeout=30):
    """Fetch a URL and return the raw bytes. Raises urllib.error.URLError
    on failure. Sets a polite User-Agent.
    """
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "SJWP-2026-NemA-pipeline/1.0 (academic use; "
                          "see LICENSE_NOTIFICATION.txt)"
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def fetch_8bpq(data_dir):
    """Fetch 8BPQ mmCIF (decompressed) and chain A as .pdb.

    Returns the path to the .cif file and the .pdb file.
    """
    os.makedirs(data_dir, exist_ok=True)
    cif_path = os.path.join(data_dir, f"{PDB_ID}.cif")
    pdb_path = os.path.join(data_dir, f"{PDB_ID}_chainA.pdb")

    if os.path.exists(cif_path):
        print(f"[fetch] {cif_path} already exists; skipping download.")
    else:
        print(f"[fetch] Downloading {PDB_ID}.cif from RCSB ...")
        try:
            raw = fetch_url(RCSB_DOWNLOAD_URL)
        except Exception as e1:
            print(f"[fetch] Direct .cif download failed: {e1}")
            print(f"[fetch] Trying .cif.gz ...")
            try:
                raw_gz = fetch_url(RCSB_DOWNLOAD_URL_GZ)
                raw = gzip.decompress(raw_gz)
            except Exception as e2:
                print(f"[fetch] .cif.gz also failed: {e2}")
                print(f"[fetch] Trying PDBj mirror ...")
                raw = fetch_url(PDBJ_FALLBACK_URL)
        with open(cif_path, "wb") as f:
            f.write(raw)
        print(f"[fetch] Wrote {cif_path} ({len(raw):,} bytes)")

    # Parse the mmCIF and extract chain A as a .pdb file
    if not os.path.exists(pdb_path):
        print(f"[parse] Extracting chain A as .pdb ...")
        try:
            from Bio.PDB import MMCIFParser, PDBIO, Select
        except ImportError:
            print(f"[parse] biopython not available; skipping chain A extraction.")
            print(f"[parse] mmCIF will be used directly in the analysis step.")
            return cif_path, None

        parser = MMCIFParser(QUIET=True)
        # biopython 1.87's MMCIFParser takes a file path, not a BytesIO
        structure = parser.get_structure(PDB_ID, cif_path)
        # Pick the first model and chain A
        model = structure[0]

        # If chain A is not present, list the chains and pick the first
        chains = [c.id for c in model]
        if "A" in chains:
            chain_id = "A"
        else:
            chain_id = sorted(chains)[0]
            print(f"[parse] Chain A not present; using chain {chain_id} instead.")

        class ChainSelect(Select):
            def accept_chain(self, chain):
                return chain.id == chain_id

        io_writer = PDBIO()
        io_writer.set_structure(structure)
        io_writer.save(pdb_path, ChainSelect())
        print(f"[parse] Wrote {pdb_path} (chain {chain_id}, {len(model[chain_id])} residues)")

    return cif_path, pdb_path


if __name__ == "__main__":
    data_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "data"
    )
    cif, pdb = fetch_8bpq(data_dir)
    print(f"[done] cif: {cif}")
    print(f"[done] pdb: {pdb}")
