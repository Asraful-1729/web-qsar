"""
Phase 1, BUILD_PLAN.md item 7: full-population recompute of the
mechanistically-supported tier (Rev 5 Section 4.1's two-of-three rule),
extending Phase 0b's fetch_mechanism_support.py (which covered only the
710 sampled query drugs) to the FULL matched population (~4,715 drugs
across Sets A+B).

Efficiency fix over the Phase 0b version: batches molecule_chembl_id__in
queries (25 molecules/request, confirmed working) instead of one request
per molecule -- the earlier per-molecule pull (710 molecules) took ~95
minutes; this should scale much better across ~4,700.

Query also filters standard_type__in=IC50,Ki,Kd,EC50 and
pchembl_value__isnull=false server-side, matching the exact ground-truth
recipe used everywhere else in this program (PHASE1_BENCHMARK_CARD.md
Section 2) -- an unfiltered first attempt pulled every activity record
regardless of type/relation (e.g. 4,107 rows for one promiscuous molecule,
CHEMBL521, vs. 134 filtered), which was both off-recipe and ~30x slower;
projected completion at that rate was ~36 hours for the full population.

Usage: python3 fetch_mechanism_support_full.py
"""
import json
import os
import sys
import time
import urllib.request
import urllib.parse

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
BATCH_SIZE = 25


def _get(url, retries=5):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                return json.load(resp)
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def fetch_activities_for_batch(mcids):
    ids_param = ",".join(mcids)
    url = (f"{BASE}/activity.json?molecule_chembl_id__in={urllib.parse.quote(ids_param)}"
           f"&standard_type__in=IC50,Ki,Kd,EC50&pchembl_value__isnull=false&limit=1000"
           f"&only=molecule_chembl_id,target_chembl_id,document_chembl_id,assay_type,target_organism,pchembl_value")
    rows = []
    while url:
        d = _get(url)
        rows.extend(d["activities"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def main():
    phase0_dir = os.path.join(HERE, "..", "phase0")
    with open(os.path.join(phase0_dir, "data", "eval_sets.json")) as f:
        eval_sets = json.load(f)

    by_mcid = {}
    for key in ["set_A_approved", "set_B_clinical"]:
        for r in eval_sets[key]:
            by_mcid[r["molecule_chembl_id"]] = r

    mcids = list(by_mcid.keys())
    print(f"{len(mcids)} distinct drugs in the full matched population (Sets A+B)", flush=True)

    out = {}
    out_path = os.path.join(DATA_DIR, "mechanism_support_full.json")
    if os.path.exists(out_path):
        with open(out_path) as f:
            out = json.load(f)
        print(f"  resuming: {len(out)} already done", flush=True)

    remaining = [m for m in mcids if m not in out]
    print(f"  {len(remaining)} remaining", flush=True)

    t0 = time.time()
    for bstart in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[bstart:bstart + BATCH_SIZE]
        try:
            activities = fetch_activities_for_batch(batch)
        except Exception as e:
            print(f"  FAILED batch starting {batch[0]}: {e}", flush=True)
            continue

        per_mol = {}
        for a in activities:
            mcid = a.get("molecule_chembl_id")
            tcid = a.get("target_chembl_id")
            if mcid is None or tcid is None:
                continue
            row = by_mcid.get(mcid)
            if row is None or tcid not in set(row["annotated_targets"]):
                continue
            entry = per_mol.setdefault(mcid, {}).setdefault(tcid, {"documents": set(), "has_binding_assay": False, "organisms": set()})
            doc = a.get("document_chembl_id")
            if doc:
                entry["documents"].add(doc)
            if a.get("assay_type") == "B":
                entry["has_binding_assay"] = True
            org = a.get("target_organism")
            if org:
                entry["organisms"].add(org)

        for mcid in batch:
            targets = per_mol.get(mcid, {})
            out[mcid] = {
                tcid: {"n_documents": len(e["documents"]), "has_binding_assay": e["has_binding_assay"],
                       "organisms_seen": sorted(e["organisms"])}
                for tcid, e in targets.items()
            }

        done = bstart + len(batch)
        elapsed = time.time() - t0
        rate = done / elapsed if elapsed > 0 else 0
        remaining_time = (len(remaining) - done) / rate if rate > 0 else 0
        print(f"  {done}/{len(remaining)} ({elapsed:.0f}s elapsed, ~{remaining_time/60:.1f} min remaining)", flush=True)

        if done % (BATCH_SIZE * 4) == 0 or done == len(remaining):
            with open(out_path, "w") as f:
                json.dump(out, f)

    with open(out_path, "w") as f:
        json.dump(out, f)
    print(f"Done: {len(out)} molecules -> {out_path} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
