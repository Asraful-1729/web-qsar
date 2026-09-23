"""
Phase 2 rebuild follow-up: document_chembl_id -> publication year backfill,
unblocking the temporal holdout that PHASE7_RELEASE_READINESS.md flagged
as the #1 remaining G1 blocker. Confirmed live: document.json exposes a
real `year` field (publication year) per document -- activity.json itself
does NOT serialize it (same class of gap as confidence_score, item 2's
backfill), so this needs the same batched separate-endpoint lookup.

54,717 distinct document_chembl_id values across Stages 1-5b's raw data
(much smaller than confidence_score's 270,225 distinct assay ids -- many
assays share a document). Same proven pattern: batched (25 ids/request via
document_chembl_id__in), 30 concurrent workers, checkpointed/resumable.

This gives YEAR-level granularity, not exact date -- sufficient for
rev5's own "renewable temporal holdout" framing (regenerable per ChEMBL
release, coarse cutoffs), not claimed as finer than it is.

Output: data/document_year_map.json -- {document_chembl_id: year_or_null}

Usage: python3 stage_document_year_backfill.py
"""
import json
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


def fetch_batch(doc_ids):
    ids_param = ",".join(doc_ids)
    url = (f"{BASE}/document.json?document_chembl_id__in={urllib.parse.quote(ids_param)}"
           f"&only=document_chembl_id,year&limit=1000")
    rows = []
    while url:
        d = _get(url)
        rows.extend(d["documents"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def collect_document_ids():
    ids = set()
    files = ["stage1_activities.jsonl", "stage2_censored.jsonl", "stage3_potency.jsonl",
              "stage4_inactives.jsonl", "stage5b_orthologue_activities.jsonl"]
    for fn in files:
        with open(os.path.join(DATA_DIR, fn)) as f:
            for line in f:
                r = json.loads(line)
                d = r.get("document_chembl_id")
                if d:
                    ids.add(d)
    return sorted(ids)


def main():
    print("Collecting distinct document_chembl_id values...", flush=True)
    all_ids = collect_document_ids()
    print(f"{len(all_ids)} distinct document ids", flush=True)

    out_path = os.path.join(DATA_DIR, "document_year_map.json")
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
                    result[r["document_chembl_id"]] = r.get("year")
                n_done_batches += 1
                if n_done_batches % 20 == 0:
                    with open(out_path, "w") as f:
                        json.dump(result, f)
            elapsed = time.time() - t0
            rate = n_done_batches * BATCH_SIZE / elapsed if elapsed > 0 else 0
            print(f"  {len(result)}/{len(all_ids)} document ids resolved "
                  f"({elapsed:.0f}s elapsed, ~{rate:.0f} ids/s)", flush=True)

    with open(out_path, "w") as f:
        json.dump(result, f)
    n_with_year = sum(1 for v in result.values() if v is not None)
    print(f"\nDone: {len(result)}/{len(all_ids)} resolved, {n_with_year} have a real year -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
