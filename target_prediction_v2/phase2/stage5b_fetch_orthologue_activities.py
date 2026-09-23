"""
Phase 2 rebuild, Stage 5b (PHASE2_REBUILD_EXECUTION_PLAN.md): fetch
activity data for exactly the 513 confirmed orthologue target ids from
Stage 5a's gene-symbol map, batched via target_chembl_id__in (25
ids/request) -- NOT a per-species unrestricted bulk pull.

REVISION from the original per-species-bulk-pull design: measured live
that an unrestricted Mus musculus pull returns 90,341 records for just 4
relevant targets (target_chembl_id__in restricted to those 4 returns the
real 7 matching records) -- pulling everything and filtering locally would
waste >99.9% of the transfer. target_chembl_id__in is the right tool here
(same pattern already proven in Phase 1's mechanism-support fetch), not
the per-species bulk approach used for the earlier stages (which made
sense there because the FILTER itself was already selective; here the
selectivity comes from a small, known target id list instead).

Same bioactivity filter as Stage 1: IC50/Ki/Kd/EC50, confidence>=8,
pchembl_value present (implies '=' relation, per L4).

Tags every row with species_provenance=<organism> and the human
target_chembl_id(s) it maps to -- kept as a separate stratum (item 4).

Usage: python3 stage5b_fetch_orthologue_activities.py
"""
import json
import os
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
BATCH_SIZE = 25

FIELDS = ("molecule_chembl_id,target_chembl_id,standard_type,standard_relation,"
          "standard_value,standard_units,pchembl_value,confidence_score,"
          "assay_type,assay_chembl_id,document_chembl_id,canonical_smiles")


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


def fetch_batch(target_ids):
    ids_param = ",".join(target_ids)
    url = (f"{BASE}/activity.json?target_chembl_id__in={urllib.parse.quote(ids_param)}"
           f"&standard_type__in=IC50,Ki,Kd,EC50&confidence_score__gte=8"
           f"&pchembl_value__isnull=false&limit=1000&only={FIELDS}")
    rows = []
    while url:
        d = _get(url)
        rows.extend(d["activities"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def main():
    with open(os.path.join(DATA_DIR, "orthologue_target_map.json")) as f:
        m = json.load(f)

    by_orthologue_target = {}
    organism_of = {}
    for human_tid, v in m.items():
        for o in v["orthologues"]:
            by_orthologue_target.setdefault(o["target_chembl_id"], []).append(human_tid)
            organism_of[o["target_chembl_id"]] = o["organism"]

    all_target_ids = list(by_orthologue_target.keys())
    print(f"{len(all_target_ids)} distinct orthologue target ids to fetch", flush=True)

    out_path = os.path.join(DATA_DIR, "stage5b_orthologue_activities.jsonl")
    ckpt_path = os.path.join(DATA_DIR, "stage5b_done_targets.json")

    done = set()
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            done = set(json.load(f))
        print(f"Resuming: {len(done)} target ids already fetched", flush=True)

    remaining = [t for t in all_target_ids if t not in done]
    t0 = time.time()
    n_rows_written = 0

    with open(out_path, "a") as f:
        for bstart in range(0, len(remaining), BATCH_SIZE):
            batch = remaining[bstart:bstart + BATCH_SIZE]
            try:
                rows = fetch_batch(batch)
            except Exception as e:
                print(f"  FAILED batch starting {batch[0]}: {e} -- will retry on next resume", flush=True)
                continue

            for r in rows:
                tid = r.get("target_chembl_id")
                if tid not in by_orthologue_target:
                    continue
                r["species_provenance"] = organism_of[tid]
                r["human_target_chembl_ids"] = by_orthologue_target[tid]
                f.write(json.dumps(r) + "\n")
                n_rows_written += 1
            f.flush()

            done.update(batch)
            with open(ckpt_path, "w") as cf:
                json.dump(sorted(done), cf)

            elapsed = time.time() - t0
            done_this_run = bstart + len(batch)
            print(f"  {len(done)}/{len(all_target_ids)} target ids done "
                  f"({elapsed:.0f}s elapsed, {n_rows_written} activity rows written)", flush=True)

    print(f"\nDone: {len(done)}/{len(all_target_ids)} target ids -> {out_path} "
          f"({n_rows_written} total rows)", flush=True)


if __name__ == "__main__":
    main()
