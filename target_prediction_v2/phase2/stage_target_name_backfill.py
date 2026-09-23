"""
Reproducibility/serving-layer follow-up (research-handoff roadmap items 4/5):
`build_v2_index_files.py` deliberately wrote target_pref_name/target_id as
None (disclosed in its own docstring: not needed for the density-adaptive
refit). It IS needed now that a live serving module is being built --
showing a caller a bare ChEMBL target id with no name is not a usable
research-tool result.

Same proven pattern as stage_document_year_backfill.py / the confidence
backfill: batched (target_chembl_id__in), 30 concurrent workers,
checkpointed/resumable. 4,882 distinct target_chembl ids in v2_index
(single-protein human targets only, per Phase 2's own verified allowlist).

Output: data/target_name_map.json -- {target_chembl_id: pref_name_or_null}

Usage: python3 stage_target_name_backfill.py
"""
import json
import os
import threading
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

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


def fetch_batch(target_ids):
    ids_param = ",".join(target_ids)
    url = (f"{BASE}/target.json?target_chembl_id__in={urllib.parse.quote(ids_param)}"
           f"&only=target_chembl_id,pref_name&limit=1000")
    rows = []
    while url:
        d = _get(url)
        rows.extend(d["targets"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def collect_target_ids():
    df = pd.read_csv(os.path.join(HERE, "v2_index", "compounds.csv.gz"), usecols=["target_chembl"])
    ortho_path = os.path.join(HERE, "..", "phase3", "orthologue_index", "compounds.csv.gz")
    ids = set(df["target_chembl"].unique())
    if os.path.exists(ortho_path):
        odf = pd.read_csv(ortho_path, usecols=["target_chembl"])
        ids |= set(odf["target_chembl"].unique())
    return sorted(ids)


def main():
    print("Collecting distinct target_chembl ids from v2_index + orthologue_index...", flush=True)
    all_ids = collect_target_ids()
    print(f"{len(all_ids)} distinct target ids", flush=True)

    out_path = os.path.join(DATA_DIR, "target_name_map.json")
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
                seen = set()
                for r in rows:
                    result[r["target_chembl_id"]] = r.get("pref_name")
                    seen.add(r["target_chembl_id"])
                for tid in batch:
                    if tid not in seen:
                        result[tid] = None
                n_done_batches += 1
                if n_done_batches % 20 == 0:
                    with open(out_path, "w") as f:
                        json.dump(result, f)
            elapsed = time.time() - t0
            rate = n_done_batches * BATCH_SIZE / elapsed if elapsed > 0 else 0
            print(f"  {len(result)}/{len(all_ids)} target ids resolved "
                  f"({elapsed:.0f}s elapsed, ~{rate:.0f} ids/s)", flush=True)

    with open(out_path, "w") as f:
        json.dump(result, f)
    n_with_name = sum(1 for v in result.values() if v)
    print(f"\nDone: {len(result)}/{len(all_ids)} resolved, {n_with_name} have a real name -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
