"""
Phase 1 item 6 -- PILOT of an AI/database-assisted evidence-gathering pass
over adjudication_candidates.csv, on a small subset (first N rows) to check
output quality before committing to all 153 rows.

Design, deliberately NOT free-text literature guessing: the two named
categories rev5 SS4.4 cares about are checked MECHANICALLY against ChEMBL's
own structured data wherever possible:

  - salt_form_gap: query ChEMBL for sibling molecule IDs sharing the same
    molecule_hierarchy.parent_chembl_id as the query drug, then check
    directly whether any sibling has an annotated (pchembl-present) activity
    against the predicted target. This is a database FACT, not an inference.

  - orthologue_or_homologue_gap: find other ChEMBL SINGLE PROTEIN targets
    with the identical pref_name (same protein, different organism -- ChEMBL
    names orthologous targets identically across species in the large
    majority of cases) and check directly whether the SAME query drug has an
    annotated activity against that non-human target. Also a database FACT.

  - Only when neither mechanical check fires: a PubMed co-occurrence search
    (drug pref_name + target pref_name/gene) as a weak, clearly-labeled
    secondary signal -- an evidence COUNT, never a verdict.

DrugBank is not used: its public pages return 403 (bot-blocked) and no API
key is available -- confirmed by a live connectivity test, not assumed.

This produces ai_* columns only. It never writes to orthologue_or_homologue_
annotation_gap / salt_form_annotation_gap / final_verdict -- those stay for
a human reviewer, per PHASE1_ADJUDICATION_STUDY_STATUS.md.

Usage: python3 adjudication_pilot.py [n_rows]
"""
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HERE = os.path.dirname(os.path.abspath(__file__))


def get_json(url, retries=4):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=20) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt == retries - 1:
                raise
            time.sleep(1.5 ** attempt)
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 ** attempt)


def molecule_info(mcid):
    m = get_json(f"{CHEMBL}/molecule/{mcid}.json")
    if not m:
        return None, None
    parent = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or mcid
    return m.get("pref_name"), parent


def sibling_salt_forms(parent_chembl_id, exclude_mcid):
    url = f"{CHEMBL}/molecule.json?molecule_hierarchy__parent_chembl_id={parent_chembl_id}&limit=50"
    d = get_json(url)
    if not d:
        return []
    return [m["molecule_chembl_id"] for m in d.get("molecules", []) if m["molecule_chembl_id"] != exclude_mcid]


def has_activity(mcid, tcid):
    url = (f"{CHEMBL}/activity.json?molecule_chembl_id={mcid}&target_chembl_id={tcid}"
           f"&pchembl_value__isnull=false&limit=1")
    d = get_json(url)
    if not d:
        return False, None
    acts = d.get("activities", [])
    if acts:
        return True, acts[0].get("pchembl_value")
    return False, None


def target_info(tcid):
    t = get_json(f"{CHEMBL}/target/{tcid}.json")
    if not t:
        return None, None, None
    gene_symbol = None
    for tc in t.get("target_components", []):
        for syn in tc.get("target_component_synonyms", []):
            if syn.get("syn_type") == "GENE_SYMBOL":
                gene_symbol = syn.get("component_synonym")
                break
        if gene_symbol:
            break
    return t.get("pref_name"), t.get("target_type"), gene_symbol


def orthologue_targets(gene_symbol, pref_name, exclude_tcid):
    """Prefer gene-symbol synonym matching (catches orthologues ChEMBL names
    identically or near-identically across species, e.g. rat/mouse FKBP1A)
    over exact pref_name matching, which misses any species-specific naming
    variation. Falls back to pref_name if no gene symbol is available."""
    if gene_symbol:
        url = f"{CHEMBL}/target.json?target_synonym__iexact={urllib.parse.quote(gene_symbol)}&limit=20"
    elif pref_name:
        url = (f"{CHEMBL}/target.json?pref_name__iexact={urllib.parse.quote(pref_name)}"
               f"&target_type={urllib.parse.quote('SINGLE PROTEIN')}&limit=20")
    else:
        return []
    d = get_json(url)
    if not d:
        return []
    return [(t["target_chembl_id"], (t.get("organism") or "?"))
            for t in d.get("targets", [])
            if t["target_chembl_id"] != exclude_tcid and t.get("target_type") == "SINGLE PROTEIN"]


def pubmed_cooccurrence(drug_name, target_name, gene_symbol=None):
    """Gene-symbol/synonym-aware fallback -- matches rev5's request to not
    rely solely on exact full-protein-name phrase matching, which is overly
    conservative (e.g. 'Kappa-type opioid receptor' rarely appears verbatim
    even in papers clearly about that receptor, while the gene symbol or a
    short common name does)."""
    if not drug_name:
        return 0, [], None
    target_terms = []
    if gene_symbol:
        target_terms.append(f'"{gene_symbol}"[Title/Abstract]')
    if target_name:
        target_terms.append(f'"{target_name}"[Title/Abstract]')
    if not target_terms:
        return 0, [], None
    target_clause = "(" + " OR ".join(target_terms) + ")"
    term = urllib.parse.quote(f'("{drug_name}"[Title/Abstract]) AND {target_clause}')
    url = f"{EUTILS}/esearch.fcgi?db=pubmed&term={term}&retmode=json&retmax=3"
    d = get_json(url)
    if not d:
        return 0, [], None
    res = d.get("esearchresult", {})
    count = int(res.get("count", 0))
    return count, res.get("idlist", []), (gene_symbol or target_name)


def evaluate_row(row):
    mcid, tcid = row["molecule_chembl_id"], row["predicted_target_chembl_id"]
    drug_name, parent = molecule_info(mcid)
    target_name, target_type, gene_symbol = target_info(tcid)

    out = {
        "ai_drug_name": drug_name,
        "ai_target_name": target_name,
        "ai_target_gene_symbol": gene_symbol,
        "ai_salt_form_gap": "No",
        "ai_salt_form_evidence": "",
        "ai_orthologue_gap": "No",
        "ai_orthologue_evidence": "",
        "ai_pubmed_cooccurrence_count": "",
        "ai_pubmed_example_pmids": "",
        "ai_pubmed_term_used": "",
        "ai_suggested_verdict": "genuine_miss (no mechanical evidence found)",
        "ai_confidence": "Low",
    }

    siblings = sibling_salt_forms(parent, mcid) if parent else []
    for sib in siblings[:10]:
        found, pchembl = has_activity(sib, tcid)
        if found:
            out["ai_salt_form_gap"] = "Yes"
            out["ai_salt_form_evidence"] = f"Sibling salt-form {sib} (same parent {parent}) has annotated activity vs {tcid} (pChEMBL={pchembl})"
            out["ai_suggested_verdict"] = "real_but_unannotated (salt-form gap)"
            out["ai_confidence"] = "High (direct ChEMBL fact)"
            break

    if out["ai_salt_form_gap"] == "No" and target_type == "SINGLE PROTEIN" and (gene_symbol or target_name):
        orthologues = orthologue_targets(gene_symbol, target_name, tcid)
        for other_tcid, organism in orthologues[:10]:
            found, pchembl = has_activity(mcid, other_tcid)
            if found:
                out["ai_orthologue_gap"] = "Yes"
                out["ai_orthologue_evidence"] = f"Same drug has annotated activity vs {other_tcid} ({organism} ortholog of {target_name}, pChEMBL={pchembl})"
                out["ai_suggested_verdict"] = "real_but_unannotated (orthologue/homologue gap)"
                out["ai_confidence"] = "High (direct ChEMBL fact)"
                break

    if out["ai_salt_form_gap"] == "No" and out["ai_orthologue_gap"] == "No":
        count, pmids, term_used = pubmed_cooccurrence(drug_name, target_name, gene_symbol)
        out["ai_pubmed_cooccurrence_count"] = count
        out["ai_pubmed_example_pmids"] = ";".join(pmids)
        out["ai_pubmed_term_used"] = term_used
        if count >= 1:
            out["ai_suggested_verdict"] = "uncertain (literature co-occurrence found, not verified)"
            out["ai_confidence"] = f"Low-Medium ({count} PubMed hits, unread)"

    return out


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    in_path = os.path.join(HERE, "data", "adjudication_candidates.csv")
    with open(in_path) as f:
        rows = list(csv.DictReader(f))
    subset = rows[:n]

    ai_fields = ["ai_drug_name", "ai_target_name", "ai_target_gene_symbol",
                 "ai_salt_form_gap", "ai_salt_form_evidence",
                 "ai_orthologue_gap", "ai_orthologue_evidence",
                 "ai_pubmed_cooccurrence_count", "ai_pubmed_example_pmids", "ai_pubmed_term_used",
                 "ai_suggested_verdict", "ai_confidence"]
    fieldnames = list(subset[0].keys()) + ai_fields

    out_path = os.path.join(HERE, "data", "adjudication_pilot_output.csv")
    done_keys = set()
    if os.path.exists(out_path):
        with open(out_path) as f:
            for r in csv.DictReader(f):
                done_keys.add((r["molecule_chembl_id"], r["predicted_target_chembl_id"]))
        print(f"Resuming: {len(done_keys)} rows already done in {out_path}", flush=True)

    file_exists = os.path.exists(out_path) and len(done_keys) > 0
    with open(out_path, "a" if file_exists else "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            w.writeheader()
        for i, row in enumerate(subset):
            key = (row["molecule_chembl_id"], row["predicted_target_chembl_id"])
            if key in done_keys:
                continue
            print(f"[{i+1}/{len(subset)}] {row['molecule_chembl_id']} -> {row['predicted_target_chembl_id']}", flush=True)
            try:
                ev = evaluate_row(row)
            except Exception as e:
                print(f"    FAILED: {e} -- skipping, will retry on next resume", flush=True)
                continue
            merged = {**row, **ev}
            w.writerow(merged)
            f.flush()
            print(f"    verdict: {ev['ai_suggested_verdict']} ({ev['ai_confidence']})", flush=True)

    with open(out_path) as f:
        results = list(csv.DictReader(f))
    n_salt = sum(1 for r in results if r["ai_salt_form_gap"] == "Yes")
    n_ortho = sum(1 for r in results if r["ai_orthologue_gap"] == "Yes")
    n_pubmed_only = sum(1 for r in results if "uncertain" in r["ai_suggested_verdict"])
    n_no_evidence = len(results) - n_salt - n_ortho - n_pubmed_only
    print(f"\nOf {len(results)} total in {out_path}: {n_salt} salt-form gap (mechanical fact), "
          f"{n_ortho} orthologue gap (mechanical fact), {n_pubmed_only} uncertain (PubMed hit only, unread), "
          f"{n_no_evidence} no evidence found (likely genuine miss)")


if __name__ == "__main__":
    main()
