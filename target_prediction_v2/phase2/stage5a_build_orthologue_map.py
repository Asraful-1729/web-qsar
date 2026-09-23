"""
Phase 2 rebuild, Stage 5a (PHASE2_REBUILD_EXECUTION_PLAN.md): build the
gene-symbol-based orthologue map -- item 4's decided method (confirmed
precise in the adjudication pilot), applied at full scale.

Design improvement over the plan's original per-species-pull-then-
intersect description: ChEMBL's target list endpoint already returns full
target_components/synonyms in the SAME bulk paginated response used to
build single_protein_target_ids.json (confirmed by inspection -- no
per-target detail calls needed). And the full cross-organism single-
protein target universe is only 11,055 targets (11 pages) -- small enough
to pull entirely in one shot (target_type=SINGLE PROTEIN, no organism
filter) and group by gene symbol locally, rather than querying human
targets and per-species orthologues separately.

Output: data/orthologue_target_map.json --
  {human_target_chembl_id: {"gene_symbol": str,
                             "orthologues": [{"target_chembl_id": str, "organism": str}, ...]}}
Only human targets with >=1 confirmed non-human orthologue (same gene
symbol, restricted to P0b-19's species list) are included.

Usage: python3 stage5a_build_orthologue_map.py
"""
import json
import math
import os
import time
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")

# P0b-19's species list (PHASE0B_ADDENDUM.md Section 10)
ORTHOLOGUE_SPECIES = {
    "Mus musculus", "Rattus norvegicus", "Bos taurus", "Sus scrofa",
    "Oryctolagus cuniculus", "Canis familiaris", "Canis lupus familiaris",
    "Cavia porcellus", "Macaca mulatta", "Gallus gallus",
}


def _get(url, retries=9):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                return json.load(resp)
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(min(2 ** attempt, 60))


def gene_symbol_of(target):
    for tc in target.get("target_components", []):
        for syn in tc.get("target_component_synonyms", []):
            if syn.get("syn_type") == "GENE_SYMBOL":
                return syn.get("component_synonym")
    return None


def main():
    url = f"{BASE}/target.json?target_type=SINGLE+PROTEIN&limit=1000"
    targets = []
    t0 = time.time()
    while url:
        d = _get(url)
        targets.extend(d["targets"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
        print(f"  {len(targets)} single-protein targets fetched so far ({time.time()-t0:.0f}s)", flush=True)

    print(f"Total: {len(targets)} single-protein targets, any organism", flush=True)

    by_gene = {}
    for t in targets:
        gs = gene_symbol_of(t)
        if not gs:
            continue
        by_gene.setdefault(gs, []).append({
            "target_chembl_id": t["target_chembl_id"],
            "organism": t.get("organism"),
        })

    out = {}
    for gene_symbol, members in by_gene.items():
        human = [m for m in members if m["organism"] == "Homo sapiens"]
        orthologues = [m for m in members if m["organism"] in ORTHOLOGUE_SPECIES]
        if human and orthologues:
            for h in human:
                out[h["target_chembl_id"]] = {
                    "gene_symbol": gene_symbol,
                    "orthologues": orthologues,
                }

    out_path = os.path.join(DATA_DIR, "orthologue_target_map.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    n_orthologue_targets = len(set(o["target_chembl_id"] for v in out.values() for o in v["orthologues"]))
    print(f"\n{len(out)} human targets have >=1 confirmed orthologue "
          f"(P0b-19 counted 2,600 candidates via a weaker name-matched proxy -- "
          f"this gene-symbol-based count is the precise figure)", flush=True)
    print(f"{n_orthologue_targets} distinct non-human orthologue target ids across "
          f"{len(ORTHOLOGUE_SPECIES)} candidate species", flush=True)
    print(f"Wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
