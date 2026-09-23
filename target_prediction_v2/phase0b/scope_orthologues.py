"""
Phase 0b, P0b-19 -- orthologue scoping (rev5's actual, narrower scope: "a
cheap counting exercise... informs Phase 2", NOT a full second evaluation
arm -- see the divergence from the pasted directive noted at the start of
this session).

Pulls ChEMBL's full single-protein target list across ALL organisms
(11,055 total, checked interactively: cheap, ~11 paginated calls), and
locally joins against v1's own 4,658 human single-protein targets by exact
`pref_name` match -- a reasonable, disclosed proxy for "orthologue" (ChEMBL
gives shared protein names to orthologous single-protein targets across
species in the overwhelming majority of cases; this is NOT the same rigor
as a HomoloGene-style formal orthology call, and is reported as a proxy,
not a validated mapping -- a real orthologue tier (Phase 2) would need the
proper mapping [23]).

Output: data/orthologue_scoping.json
"""
import json
import os
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE0_DIR = os.path.join(HERE, "..", "phase0")
BASE = "https://www.ebi.ac.uk/chembl/api/data"


def _get(url, retries=5):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return json.load(resp)
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def fetch_all_single_protein_targets():
    rows = []
    url = f"{BASE}/target.json?target_type=SINGLE%20PROTEIN&limit=1000&only=target_chembl_id,pref_name,organism"
    while url:
        d = _get(url)
        rows.extend(d["targets"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def main():
    import sys
    sys.path.insert(0, PHASE0_DIR)
    sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
    os.environ.setdefault("TARGET_FISHING_INDEX_DIR", os.path.abspath(
        os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "target_fishing_index")))
    import target_fishing as TF
    _, _, full_df = TF._load()
    human_targets = full_df.drop_duplicates(subset="target_chembl")[["target_chembl", "target_pref_name"]]
    human_pref_names = set(n for n in human_targets["target_pref_name"].dropna().tolist())
    print(f"v1 index: {len(human_targets)} distinct human targets, {len(human_pref_names)} distinct pref_names", flush=True)

    print("Fetching all single-protein targets, all organisms...", flush=True)
    t0 = time.time()
    all_targets = fetch_all_single_protein_targets()
    print(f"  {len(all_targets)} total ({time.time()-t0:.0f}s)", flush=True)

    by_organism = {}
    orthologue_candidates = []
    for t in all_targets:
        org = t.get("organism") or "unknown"
        by_organism[org] = by_organism.get(org, 0) + 1
        if org != "Homo sapiens" and t.get("pref_name") in human_pref_names:
            orthologue_candidates.append({"target_chembl_id": t["target_chembl_id"], "pref_name": t["pref_name"], "organism": org})

    human_pref_names_covered = set(c["pref_name"] for c in orthologue_candidates)
    org_counts_candidates = {}
    for c in orthologue_candidates:
        org_counts_candidates[c["organism"]] = org_counts_candidates.get(c["organism"], 0) + 1

    out = {
        "n_human_single_protein_targets_v1": len(human_targets),
        "n_distinct_human_pref_names_v1": len(human_pref_names),
        "n_all_single_protein_targets_all_organisms": len(all_targets),
        "n_organisms_represented": len(by_organism),
        "top_20_organisms_by_target_count": sorted(by_organism.items(), key=lambda kv: -kv[1])[:20],
        "n_orthologue_candidate_targets_by_pref_name_match": len(orthologue_candidates),
        "n_distinct_human_pref_names_with_a_orthologue_candidate": len(human_pref_names_covered),
        "pct_human_targets_with_a_orthologue_candidate": round(100 * len(human_pref_names_covered) / len(human_pref_names), 1),
        "orthologue_candidate_organism_breakdown": sorted(org_counts_candidates.items(), key=lambda kv: -kv[1])[:20],
        "method_caveat": "pref_name exact-match proxy for orthology, not a validated HomoloGene-style mapping; "
                          "does not yet count how many ADDITIONAL compounds/annotations each candidate target would "
                          "bring (that needs a per-target activity pull, out of scope for this counting pass).",
    }
    out_path = os.path.join(HERE, "data", "orthologue_scoping.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2)[:2000], flush=True)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
