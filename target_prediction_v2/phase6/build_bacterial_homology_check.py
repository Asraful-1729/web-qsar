"""
Phase 6 scoping: human-homology annotation feasibility, reusing the exact
gene-symbol-matching method already validated for the orthologue tier
(phase2/stage5a_build_orthologue_map.py, phase3's adjudication pilot).

Fetches full target records (with target_components/gene symbols) for the
665 bacterial single-protein targets AND the 5,869 human single-protein
targets, then checks: does a bacterial target's gene symbol also appear
among human single-protein targets? A match is a real, disclosed PROXY
for "has a close human homolog" (same caveat as the orthologue tier's
proxy nature) -- exact gene-symbol identity across kingdoms is rare and
would typically indicate a genuinely conserved/essential gene family
(e.g. DNA gyrase, ribosomal proteins), which is directly relevant to
antibacterial target selection (a target good match increases off-target
toxicity risk; a clean non-match is favorable).

Usage: python3 build_bacterial_homology_check.py
"""
import json
import os
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
PHASE2_DATA = os.path.join(HERE, "..", "phase2", "data")
BASE = "https://www.ebi.ac.uk/chembl/api/data"
BATCH_SIZE = 25


def _get(url, retries=9):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(min(2 ** attempt, 30))


def gene_symbol_of(target):
    for tc in target.get("target_components", []):
        for syn in tc.get("target_component_synonyms", []):
            if syn.get("syn_type") == "GENE_SYMBOL":
                return syn.get("component_synonym")
    return None


def fetch_gene_symbols(target_ids, cache_path):
    result = {}
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            result = json.load(f)
        print(f"  resuming: {len(result)} already fetched", flush=True)
    remaining = [t for t in target_ids if t not in result]
    t0 = time.time()
    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i:i + BATCH_SIZE]
        ids_param = urllib.parse.quote(",".join(batch))
        url = f"{BASE}/target.json?target_chembl_id__in={ids_param}&limit=1000"
        d = _get(url)
        for t in d["targets"]:
            result[t["target_chembl_id"]] = gene_symbol_of(t)
        for t in batch:
            result.setdefault(t, None)
        if (i + BATCH_SIZE) % 500 == 0 or i + BATCH_SIZE >= len(remaining):
            with open(cache_path, "w") as f:
                json.dump(result, f)
            print(f"  {min(i+BATCH_SIZE, len(remaining))}/{len(remaining)} fetched this run "
                  f"({time.time()-t0:.0f}s)", flush=True)
    with open(cache_path, "w") as f:
        json.dump(result, f)
    return result


def main():
    bacterial_targets = json.load(open(os.path.join(DATA_DIR, "bacterial_single_protein_targets.json")))
    bacterial_ids = [t["target_chembl_id"] for t in bacterial_targets]
    print(f"{len(bacterial_ids)} bacterial targets", flush=True)

    human_ids = json.load(open(os.path.join(PHASE2_DATA, "single_protein_target_ids.json")))
    print(f"{len(human_ids)} human targets", flush=True)

    print("Fetching bacterial target gene symbols...", flush=True)
    bacterial_genes = fetch_gene_symbols(bacterial_ids, os.path.join(DATA_DIR, "bacterial_gene_symbols.json"))

    print("Fetching human target gene symbols (largest pull, ~235 batches)...", flush=True)
    human_genes = fetch_gene_symbols(human_ids, os.path.join(DATA_DIR, "human_gene_symbols.json"))

    human_gene_set = {g for g in human_genes.values() if g}
    print(f"{len(human_gene_set)} distinct human gene symbols", flush=True)

    matches = []
    for tid, gene in bacterial_genes.items():
        if gene and gene.upper() in {g.upper() for g in human_gene_set}:
            matches.append({"target_chembl_id": tid, "gene_symbol": gene})

    print(f"\n{len(matches)}/{len(bacterial_ids)} bacterial targets ({100*len(matches)/len(bacterial_ids):.1f}%) "
          f"share a gene symbol with a human single-protein target "
          f"(proxy for 'has a close human homolog')", flush=True)
    print("Sample matches:", matches[:10], flush=True)

    out = {
        "n_bacterial_targets": len(bacterial_ids),
        "n_with_human_gene_symbol_match": len(matches),
        "pct_with_human_gene_symbol_match": round(100 * len(matches) / len(bacterial_ids), 1),
        "matches": matches,
    }
    out_path = os.path.join(DATA_DIR, "bacterial_human_homology_proxy.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
