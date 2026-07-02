"""structural_pipeline package — NemA*2+ structural evaluation for the
SJWP 2026 project. See individual module docstrings for details.

Modules:
  fetch_nemA_structure.py    — download PDB 8BPQ and extract chain A
  analyze_active_site.py     — FMN pocket and UniProt-annotated residues
  scan_contact_disruption.py — dimensionless contact-disruption score
  literature_check.py        — pre-docking literature sanity check
  prepare_docking.py         — prepare WT and mutant .pdbqt for Vina
  dock_chromate.py           — run Vina docking and analyze poses
  utils.py                   — shared biopython helpers

Note: this pipeline is a STRUCTURAL CONSISTENCY CHECK on the engineering
hypothesis motivating NemA*2+, not a kinetic prediction. See
Journal_Manuscript_2026.tex §3.6 and §4.3 for the full framing.
"""
