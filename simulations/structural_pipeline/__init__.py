"""structural_pipeline package — NemA*2+ structural evaluation for the
SJWP 2026 project. See individual module docstrings for details.

Modules:
  fetch_nemA_structure.py    — download PDB 8BPQ and extract chain A
  literature_check.py        — pre-docking PubMed sanity check (decision gate)
  analyze_active_site.py     — FMN pocket and UniProt-annotated residues
  prepare_docking.py         — prepare WT and ILE328ALA .pdbqt for Vina
  scan_contact_disruption.py — dimensionless contact-disruption score
  geometric_placement.py     — scipy.optimize placement of CrO4(2-)
                               (Vina fallback; Cr is not a built-in Vina type)
  utils.py                   — shared biopython helpers

Note: this pipeline is a STRUCTURAL CONSISTENCY CHECK on the engineering
hypothesis motivating NemA*2+, not a kinetic prediction. See
Journal_Manuscript_2026.tex §3.6 and §4.3 for the full framing.
"""
