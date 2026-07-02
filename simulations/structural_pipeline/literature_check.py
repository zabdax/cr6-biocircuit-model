"""
literature_check.py — Pre-docking literature sanity check for the
SJWP 2026 NemA*2+ structural pipeline.

Before running AutoDock Vina on CrO4(2-) docked into the FMN pocket
of NemA, this script queries PubMed to verify whether:

  (Q1) docking chromate (or Cr(VI)) into an FMN-containing
       flavoprotein active site has been published before
  (Q2) the mechanism of Cr(VI) reduction by NemA / OYE-family
       flavoproteins is reported as direct hydride transfer from FMN
       (which would support docking at FMN C4a) or as indirect
       (e.g., flavin-mediated ROS generation)

The decision rule:
  - If >= 2 published papers report docking chromate / Cr(VI) into
    an FMN-containing flavoprotein active site, the pipeline
    proceeds with Vina framed as a "structural compatibility"
    analysis, with Vina scores reported in the main manuscript
    text alongside geometric observables.
  - Otherwise (0-1 papers, or mechanism reported as indirect),
    the pipeline reframes the docking as an "exploratory
    structural-positioning" analysis. Vina scores are demoted to
    a supplementary CSV; only geometric observables (Cr-FMN
    C4a distance, H-bond donor counts) appear in the main text.

Output: simulations/results/literature_check.md
"""

import json
import os
import sys
import time
import urllib.parse
import urllib.request


NCBI_EMAIL = "zubayerhasanshaad99@gmail.com"  # per NCBI E-utilities policy
NCBI_TOOL = "SJWP-2026-NemA-pipeline"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _esearch(term, retmax=20):
    """Run esearch and return list of PMIDs."""
    url = (
        f"{EUTILS_BASE}/esearch.fcgi?db=pubmed&term="
        f"{urllib.parse.quote(term)}&retmax={retmax}&retmode=json"
        f"&email={NCBI_EMAIL}&tool={NCBI_TOOL}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": f"{NCBI_TOOL}/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    return data.get("esearchresult", {}).get("idlist", [])


def _esummary(pmids):
    """Run esummary for a list of PMIDs. Returns {pmid: {title, authors,
    source, pubdate, ...}}."""
    if not pmids:
        return {}
    ids = ",".join(pmids)
    url = (
        f"{EUTILS_BASE}/esummary.fcgi?db=pubmed&id={ids}&retmode=json"
        f"&email={NCBI_EMAIL}&tool={NCBI_TOOL}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": f"{NCBI_TOOL}/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    result = data.get("result", {})
    out = {}
    for pmid in pmids:
        uids = result.get(pmid, {})
        if not uids or "error" in uids:
            continue
        # Authors are stored as a list of dicts
        authors = uids.get("authors", [])
        first_author = authors[0]["name"] if authors else "?"
        out[pmid] = {
            "title": uids.get("title", ""),
            "first_author": first_author,
            "source": uids.get("source", ""),
            "pubdate": uids.get("pubdate", ""),
        }
    return out


def _efetch_abstracts(pmids):
    """Fetch abstracts as plain text for a list of PMIDs."""
    if not pmids:
        return {}
    ids = ",".join(pmids)
    url = (
        f"{EUTILS_BASE}/efetch.fcgi?db=pubmed&id={ids}&rettype=abstract"
        f"&retmode=text&email={NCBI_EMAIL}&tool={NCBI_TOOL}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": f"{NCNCBI_TOOL}/1.0" if False else f"{NCBI_TOOL}/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    # Split by PMID marker
    out = {}
    blocks = text.split("PMID:")
    for blk in blocks[1:]:
        lines = blk.strip().splitlines()
        if not lines:
            continue
        pmid_line = lines[0].strip()
        pmid = pmid_line.split()[0]
        # The abstract body follows
        body = "\n".join(lines[1:]).strip()
        out[pmid] = body
    return out


def main():
    """Run the literature queries and write the markdown report."""
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"
    )
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "literature_check.md")

    report_lines = []
    report_lines.append("# Pre-Docking Literature Sanity Check\n")
    report_lines.append("Pipeline: SJWP 2026 NemA*2+ structural evaluation")
    report_lines.append(f"Date: 2026-07-01\n")

    # ----------------------------------------------------------------
    # Q1: Docking chromate / Cr(VI) into an FMN-containing
    #     flavoprotein active site
    # ----------------------------------------------------------------
    report_lines.append("## Q1. Has docking chromate / Cr(VI) into an FMN-containing flavoprotein active site been published?\n")

    q1_queries = [
        'chromate FMN docking',
        '"chromate" AND "FMN" AND ("docking" OR "molecular docking")',
        'Cr(VI) flavoprotein docking',
        'Cr(VI) reductase docking',
        '"Old Yellow Enzyme" chromate docking',
    ]
    q1_hits = {}
    for q in q1_queries:
        try:
            pmids = _esearch(q, retmax=20)
            report_lines.append(f"- Query `{q}`: {len(pmids)} PMIDs")
            for p in pmids:
                if p not in q1_hits:
                    q1_hits[p] = q
            time.sleep(0.4)  # NCBI rate limit: 3 req/s without API key
        except Exception as e:
            report_lines.append(f"- Query `{q}`: ERROR {e}")

    report_lines.append(f"\nTotal unique PMIDs across Q1 queries: {len(q1_hits)}")
    if q1_hits:
        report_lines.append("\nTop hits (first 10):")
        meta = _esummary(list(q1_hits.keys())[:10])
        for pmid, q in list(q1_hits.items())[:10]:
            m = meta.get(pmid, {})
            title = m.get("title", "?")
            first_author = m.get("first_author", "?")
            year = m.get("pubdate", "?")[:4]
            source = m.get("source", "?")
            report_lines.append(f"- PMID {pmid} ({year}, {first_author}, {source}): {title}")

    # Inspect the abstracts of the top hits for actual docking experiments
    q1_docking_papers = 0
    q1_abstracts = {}
    if q1_hits:
        try:
            q1_abstracts = _efetch_abstracts(list(q1_hits.keys())[:5])
            time.sleep(0.4)
        except Exception as e:
            report_lines.append(f"\nWARNING: could not fetch abstracts: {e}")

    for pmid, body in q1_abstracts.items():
        body_lower = body.lower()
        # Heuristic: a docking paper will mention "dock" and one of:
        # "FMN" or "flavo" or "flavoprotein"
        if "dock" in body_lower and ("fmn" in body_lower or "flavo" in body_lower or "flavoprotein" in body_lower):
            q1_docking_papers += 1
            report_lines.append(f"  - PMID {pmid}: appears to be a docking study of chromate/Cr(VI) at a flavoprotein active site")
        elif "dock" in body_lower:
            report_lines.append(f"  - PMID {pmid}: mentions docking but not necessarily at FMN — abstract excerpt: {body[:200]!r}")
        else:
            report_lines.append(f"  - PMID {pmid}: abstract does not mention docking — likely a kinetic / spectroscopic study")

    report_lines.append(f"\n**Q1 answer: {q1_docking_papers} published paper(s) appear to report docking chromate/Cr(VI) into an FMN-containing flavoprotein active site.**")

    # ----------------------------------------------------------------
    # Q2: Mechanism of Cr(VI) reduction by NemA / OYE flavoproteins
    # ----------------------------------------------------------------
    report_lines.append("\n## Q2. Is the mechanism of Cr(VI) reduction by NemA / OYE-family flavoproteins direct FMN hydride transfer, or indirect (e.g., ROS)?\n")

    q2_queries = [
        'NemA chromate reduction mechanism',
        'chromate reductase flavoprotein mechanism',
        '"Old Yellow Enzyme" chromium reduction',
        'flavoprotein Cr(VI) hydride transfer',
    ]
    q2_hits = {}
    for q in q2_queries:
        try:
            pmids = _esearch(q, retmax=20)
            report_lines.append(f"- Query `{q}`: {len(pmids)} PMIDs")
            for p in q2_hits:
                pass
            for p in pmids:
                if p not in q2_hits:
                    q2_hits[p] = q
            time.sleep(0.4)
        except Exception as e:
            report_lines.append(f"- Query `{q}`: ERROR {e}")

    report_lines.append(f"\nTotal unique PMIDs across Q2 queries: {len(q2_hits)}")
    if q2_hits:
        report_lines.append("\nTop hits (first 10):")
        meta = _esummary(list(q2_hits.keys())[:10])
        for pmid in list(q2_hits.keys())[:10]:
            m = meta.get(pmid, {})
            title = m.get("title", "?")
            first_author = m.get("first_author", "?")
            year = m.get("pubdate", "?")[:4]
            source = m.get("source", "?")
            report_lines.append(f"- PMID {pmid} ({year}, {first_author}, {source}): {title}")

    # Inspect abstracts for the mechanism
    q2_direct_fmn = 0
    q2_indirect = 0
    q2_abstracts = {}
    if q2_hits:
        try:
            q2_abstracts = _efetch_abstracts(list(q2_hits.keys())[:5])
            time.sleep(0.4)
        except Exception as e:
            report_lines.append(f"\nWARNING: could not fetch abstracts: {e}")

    for pmid, body in q2_abstracts.items():
        body_lower = body.lower()
        if "hydride" in body_lower and "fmn" in body_lower:
            q2_direct_fmn += 1
        if "reactive oxygen" in body_lower or "ros" in body_lower or "fenton" in body_lower or "radical" in body_lower:
            q2_indirect += 1
        report_lines.append(f"  - PMID {pmid}: {body[:200]!r}")

    report_lines.append(f"\n**Q2 answer: {q2_direct_fmn} abstract(s) mention direct FMN hydride transfer; {q2_indirect} abstract(s) mention ROS / radical / Fenton-like mechanisms.**")

    # ----------------------------------------------------------------
    # Decision rule
    # ----------------------------------------------------------------
    report_lines.append("\n## Decision\n")

    if q1_docking_papers >= 2 and q2_direct_fmn >= 1 and q2_direct_fmn > q2_indirect:
        decision = "structural_compatibility"
        report_lines.append(
            "Outcome: **structural_compatibility**.\n\n"
            "Rationale: At least two published papers report docking "
            "chromate/Cr(VI) into an FMN-containing flavoprotein active site, "
            "and the prevailing mechanistic view (based on the abstracts "
            "examined) is direct FMN hydride transfer. The structural "
            "pipeline proceeds with AutoDock Vina as a structural-"
            "compatibility / substrate-positioning analysis. Vina scores "
            "are reported in the main manuscript text alongside geometric "
            "observables, with the explicit caveat that Vina's empirical "
            "scoring function is not a thermodynamic binding free energy "
            "and is reported for relative comparison only."
        )
    else:
        decision = "exploratory_positioning"
        if q1_docking_papers < 2:
            report_lines.append(
                f"- Q1 returned only {q1_docking_papers} docking precedent(s), "
                "below the >= 2 threshold for the structural-compatibility branch."
            )
        if q2_direct_fmn < 1:
            report_lines.append(
                f"- Q2 found {q2_direct_fmn} abstracts mentioning direct FMN "
                "hydride transfer; the direct-hydride mechanism is not "
                "well-supported in the surveyed abstracts."
            )
        if q2_indirect > 0:
            report_lines.append(
                f"- Q2 found {q2_indirect} abstracts mentioning ROS / Fenton-like "
                "mechanisms. If Cr(VI) reduction is not direct FMN hydride "
                "transfer, the geometric positioning of CrO4(2-) relative to "
                "FMN is not the catalytically relevant step."
            )
        report_lines.append(
            "\nOutcome: **exploratory_positioning**.\n\n"
            "Rationale: The published literature does not provide strong "
            "precedent for docking chromate into FMN-containing flavoprotein "
            "active sites, and/or the mechanism of Cr(VI) reduction is "
            "reported as not strictly direct FMN hydride transfer. The "
            "structural pipeline reframes the docking as an exploratory "
            "structural-positioning analysis. **Vina scores are demoted to "
            "a supplementary CSV** and are not used to draw conclusions in "
            "the main manuscript text. Only geometric observables (Cr-FMN "
            "C4a distance, H-bond donor counts) appear in the main results."
        )

    report_lines.append(
        f"\n## Provenance log\n"
        f"- Decision: `{decision}`\n"
        f"- Q1 docking papers found: {q1_docking_papers}\n"
        f"- Q2 direct-FMN papers: {q2_direct_fmn}\n"
        f"- Q2 ROS/Fenton papers: {q2_indirect}\n"
        f"- Q1 unique PMIDs: {len(q1_hits)}\n"
        f"- Q2 unique PMIDs: {len(q2_hits)}\n"
        f"- This file is generated by simulations/structural_pipeline/literature_check.py\n"
    )

    # Also save the decision to a machine-readable JSON
    json_path = os.path.join(out_dir, "literature_check.json")
    payload = {
        "decision": decision,
        "q1_docking_papers": q1_docking_papers,
        "q2_direct_fmn_papers": q2_direct_fmn,
        "q2_indirect_papers": q2_indirect,
        "q1_unique_pmids": len(q1_hits),
        "q2_unique_pmids": len(q2_hits),
        "q1_pmids": list(q1_hits.keys())[:20],
        "q2_pmids": list(q2_hits.keys())[:20],
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    with open(out_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")

    print(f"[literature_check] Wrote {out_path}")
    print(f"[literature_check] Wrote {json_path}")
    print(f"[literature_check] DECISION: {decision}")
    return decision, payload


if __name__ == "__main__":
    decision, payload = main()
