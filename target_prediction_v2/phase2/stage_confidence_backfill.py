"""
Phase 2 rebuild -- confidence_score backfill, found necessary after Stage 7:
ChEMBL's activity.json endpoint accepts confidence_score as a FILTER
parameter (correctly applied server-side throughout Stages 1-5b, confirmed
via exact total_count matches) but does NOT serialize confidence_score in
its own output records at all -- confirmed empirically (a full, unrestricted
activity record has no confidence_score key whatsoever). The actual value
lives on the associated Assay record instead, reachable via
assay_chembl_id (present on every activity row) -> assay.json.

270,225 distinct assay_chembl_id values across Stages 1-5b's raw data.
Batched (25 ids/request via assay_chembl_id__in, the proven Phase-1-style
pattern) with the same 30-worker concurrent + per-batch-checkpoint design
validated in fetch_common.py.

Output: data/assay_confidence_map.json -- {assay_chembl_id: confidence_score}
Stage 7's aggregation is re-run after this to properly populate the
confidence_score field it currently always writes as null.

Usage: python3 stage_confidence_backfill.py
"""
import json
import math
import os
import threading
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
BATCH_SIZE = 25
N_WORKERS = 30


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


def fetch_batch(assay_ids):
    ids_param = ",".join(assay_ids)
    url = (f"{BASE}/assay.json?assay_chembl_id__in={urllib.parse.quote(ids_param)}"
           f"&only=assay_chembl_id,confidence_score&limit=1000")
    rows = []
    while url:
        d = _get(url)
        rows.extend(d["assays"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def collect_assay_ids():
    ids = set()
    files = ["stage1_activities.jsonl", "stage2_censored.jsonl", "stage3_potency.jsonl",
              "stage4_inactives.jsonl", "stage5b_orthologue_activities.jsonl"]
    for fn in files:
        with open(os.path.join(DATA_DIR, fn)) as f:
            for line in f:
                r = json.loads(line)
                aid = r.get("assay_chembl_id")
                if aid:
                    ids.add(aid)
    return sorted(ids)


def main():
    print("Collecting distinct assay_chembl_id values...", flush=True)
    all_ids = collect_assay_ids()
    print(f"{len(all_ids)} distinct assay ids", flush=True)

    out_path = os.path.join(DATA_DIR, "assay_confidence_map.json")
    result = {}
    if os.path.exists(out_path):
        with open(out_path) as f:
            result = json.load(f)
        print(f"Resuming: {len(result)} already fetched", flush=True)

    remaining_ids = [i for i in all_ids if i not in result]
    batches = [remaining_ids[i:i + BATCH_SIZE] for i in range(0, len(remaining_ids), BATCH_SIZE)]
    print(f"{len(batches)} batches remaining", flush=True)

    write_lock = threading.Lock()
    t0 = time.time()
    n_done_batches = 0

    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(fetch_batch, b): tuple(b) for b in batches}
        for fut in as_completed(futures):
            batch = futures[fut]
            try:
                rows = fut.result()
            except Exception as e:
                print(f"  FAILED batch starting {batch[0]}: {e} -- will retry on next resume", flush=True)
                continue
            with write_lock:
                for r in rows:
                    result[r["assay_chembl_id"]] = r.get("confidence_score")
                n_done_batches += 1
                if n_done_batches % 20 == 0:
                    with open(out_path, "w") as f:
                        json.dump(result, f)
            elapsed = time.time() - t0
            rate = n_done_batches * BATCH_SIZE / elapsed if elapsed > 0 else 0
            print(f"  {len(result)}/{len(all_ids)} assay ids resolved "
                  f"({elapsed:.0f}s elapsed, ~{rate:.0f} ids/s)", flush=True)

    with open(out_path, "w") as f:
        json.dump(result, f)
    print(f"\nDone: {len(result)}/{len(all_ids)} -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
