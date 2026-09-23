"""
Phase 0b -- targeted per-record ChEMBL pull, scoped to just our sampled
query drugs (not the full 857k-compound index -- that stays aggregated,
per v1's frozen build). Fills the one real data gap flagged before Phase 0b
started: v1's index has no document_chembl_id / assay_type / standard_relation
per (compound, target) pair, so Rev 5's "mechanistically supported" tier
(Section 4.1's two-of-three rule: >=2 independent documents, a
drug_mechanism record, OR binding-assay evidence) cannot be computed from
the existing index alone.

For each sampled query drug (same seed=42 sample as capture.py), pulls
EVERY raw activity row for that molecule_chembl_id (paginated,
molecule_chembl_id is highly selective so this is cheap per molecule) and
locally aggregates, per (molecule, target) pair restricted to that drug's
own already-known annotated_targets:
  - n_distinct_documents (source publications/patents reporting this pair)
  - has_binding_assay (any row with assay_type == 'B')
  - target_organism values seen (bonus: informs P0b-19's orthologue scoping
    for free, since a non-human target_chembl_id occasionally shows up
    against the same molecule)

This does NOT re-derive annotated_targets itself (that stays exactly what
build_eval_sets.py already computed from v1's own frozen index -- changing
it here would silently redefine v1's own ground truth, which Rev 5's rules
forbid). It only adds evidence-STRENGTH metadata for targets already in
that list.

Output: data/mechanism_support.json -- {molecule_chembl_id: {target_chembl_id:
  {n_documents, has_binding_assay, organisms_seen}}}

Usage: python3 fetch_mechanism_support.py
"""
import json
import os
import sys
import time
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0"))

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")


def _get(url, retries=5):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return json.load(resp)
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def fetch_activities_for_molecule(mcid):
    rows = []
    url = (f"{BASE}/activity.json?molecule_chembl_id={mcid}&limit=1000"
           f"&only=target_chembl_id,document_chembl_id,assay_type,target_organism,pchembl_value")
    while url:
        d = _get(url)
        rows.extend(d["activities"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def main():
    import random
    with open(os.path.join(HERE, "..", "phase0", "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)

    all_query_rows = []
    for key, n, seed in [("set_A_approved", 200, 42), ("set_B_clinical", 200, 42), ("set_D_sparse_target", -1, 42)]:
        rows = list(eval_sets[key])
        rng = random.Random(seed)
        rng.shuffle(rows)
        if n > 0:
            rows = rows[:n]
        all_query_rows.extend(rows)

    by_mcid = {r["molecule_chembl_id"]: r for r in all_query_rows}
    mcids = list(by_mcid.keys())
    print(f"{len(mcids)} distinct query drugs to pull raw activity data for", flush=True)

    out = {}
    t0 = time.time()
    for i, mcid in enumerate(mcids):
        row = by_mcid[mcid]
        annotated = set(row["annotated_targets"])
        try:
            activities = fetch_activities_for_molecule(mcid)
        except Exception as e:
            print(f"  FAILED {mcid}: {e}", flush=True)
            continue

        per_target = {}
        for a in activities:
            tcid = a.get("target_chembl_id")
            if tcid not in annotated:
                continue
            entry = per_target.setdefault(tcid, {"documents": set(), "has_binding_assay": False, "organisms": set()})
            doc = a.get("document_chembl_id")
            if doc:
                entry["documents"].add(doc)
            if a.get("assay_type") == "B":
                entry["has_binding_assay"] = True
            org = a.get("target_organism")
            if org:
                entry["organisms"].add(org)

        out[mcid] = {
            tcid: {
                "n_documents": len(e["documents"]),
                "has_binding_assay": e["has_binding_assay"],
                "organisms_seen": sorted(e["organisms"]),
            }
            for tcid, e in per_target.items()
        }

        if (i + 1) % 50 == 0 or (i + 1) == len(mcids):
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            remaining = (len(mcids) - (i + 1)) / rate if rate > 0 else 0
            print(f"  {i+1}/{len(mcids)} ({elapsed:.0f}s elapsed, ~{remaining/60:.1f} min remaining)", flush=True)
            with open(os.path.join(DATA_DIR, "mechanism_support.json"), "w") as f:
                json.dump(out, f)

    with open(os.path.join(DATA_DIR, "mechanism_support.json"), "w") as f:
        json.dump(out, f)
    print(f"Done: {len(out)} molecules -> data/mechanism_support.json ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
