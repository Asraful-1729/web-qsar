"""
Phase 2 rebuild, Stage 1 (PHASE2_REBUILD_EXECUTION_PLAN.md): bulk pull of
the base '=' -relation positive evidence (same filter v1's own build used:
human, IC50/Ki/Kd/EC50, confidence>=8, pchembl_value present) -- this is
the one case where ChEMBL's own precomputed pchembl_value is safe to use
directly (per the L4 finding: it's only ever populated for '=' relations).

REVISION -- concurrent, not sequential: a first sequential pilot measured
only ~54 rec/s steady-state (implying ~13.3hr for the full ~2.6M-record
pull, ~5x worse than this plan's original v1-history-based estimate). A
live microbenchmark of concurrent requests found near-linear scaling up to
~10 workers and continued (sub-linear) gains to 30-40, with zero errors or
rate-limit responses observed at any tested concurrency -- 30 workers
measured ~585 rec/s, cutting the full pull to roughly ~70 minutes. Chose
30 workers (not 40+): captures the large majority of the achievable
speedup without needlessly maximizing concurrent load against a shared
public API.

Checkpointing is per-page-offset (not a single cursor, since pages now
complete out of order under concurrency) -- resumable from the start, the
lesson from Phase 1's mechanism-support fetch which crashed once before
checkpointing was added.

Usage:
    python3 stage1_fetch_base_evidence.py --pilot 20000   # pilot run
    python3 stage1_fetch_base_evidence.py                 # full run (resumable)
"""
import argparse
import json
import math
import os
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
PAGE_LIMIT = 1000
N_WORKERS = 30

FIELDS = ("molecule_chembl_id,target_chembl_id,standard_type,standard_relation,"
          "standard_value,standard_units,pchembl_value,confidence_score,"
          "assay_type,assay_chembl_id,document_chembl_id,canonical_smiles")

BASE_QUERY = (f"{BASE}/activity.json?target_organism=Homo+sapiens"
              f"&standard_type__in=IC50,Ki,Kd,EC50&confidence_score__gte=8"
              f"&pchembl_value__isnull=false&only={FIELDS}")


def _get(url, retries=9):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                return json.load(resp)
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(min(2 ** attempt, 60))


def fetch_page(offset):
    url = f"{BASE_QUERY}&limit={PAGE_LIMIT}&offset={offset}"
    d = _get(url)
    return offset, d["activities"]


def load_completed(ckpt_path):
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            return set(json.load(f))
    return set()


def save_completed(ckpt_path, completed):
    with open(ckpt_path, "w") as f:
        json.dump(sorted(completed), f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=int, default=None,
                     help="cap total records fetched, for a measured-throughput test run")
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    args = ap.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)
    suffix = "_pilot" if args.pilot else ""
    out_path = os.path.join(DATA_DIR, f"stage1_activities{suffix}.jsonl")
    ckpt_path = os.path.join(DATA_DIR, f"stage1_completed_offsets{suffix}.json")

    probe = _get(f"{BASE_QUERY}&limit=1")
    total_count = probe["page_meta"]["total_count"]
    print(f"Total matching records available: {total_count}", flush=True)

    n_pages = math.ceil(total_count / PAGE_LIMIT)
    all_offsets = [i * PAGE_LIMIT for i in range(n_pages)]
    if args.pilot:
        n_pilot_pages = math.ceil(args.pilot / PAGE_LIMIT)
        all_offsets = all_offsets[:n_pilot_pages]
        print(f"PILOT MODE: capping at {len(all_offsets)} pages (~{args.pilot} records)", flush=True)

    completed = load_completed(ckpt_path)
    remaining = [o for o in all_offsets if o not in completed]
    if completed:
        print(f"Resuming: {len(completed)}/{len(all_offsets)} pages already done, "
              f"{len(remaining)} remaining", flush=True)

    t0 = time.time()
    write_lock = threading.Lock()
    pages_done_this_run = 0

    with open(out_path, "a") as f:
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(fetch_page, o): o for o in remaining}
            for fut in as_completed(futures):
                offset = futures[fut]
                try:
                    off, rows = fut.result()
                except Exception as e:
                    print(f"  FAILED offset={offset}: {e} -- will retry on next resume", flush=True)
                    continue

                with write_lock:
                    for r in rows:
                        f.write(json.dumps(r) + "\n")
                    f.flush()
                    completed.add(off)
                    pages_done_this_run += 1
                    if pages_done_this_run % 10 == 0 or len(completed) == len(all_offsets):
                        save_completed(ckpt_path, completed)

                elapsed = time.time() - t0
                rec_rate = pages_done_this_run * PAGE_LIMIT / elapsed if elapsed > 0 else 0
                print(f"  {len(completed)}/{len(all_offsets)} pages "
                      f"({elapsed:.0f}s elapsed this run, ~{rec_rate:.0f} rec/s)", flush=True)

    save_completed(ckpt_path, completed)
    elapsed = time.time() - t0
    print(f"\nDone: {len(completed)}/{len(all_offsets)} pages -> {out_path} ({elapsed:.0f}s)", flush=True)
    if args.pilot and elapsed > 0 and pages_done_this_run > 0:
        measured_rate = pages_done_this_run * PAGE_LIMIT / elapsed
        est_hours = total_count / measured_rate / 3600
        print(f"Measured rate this run: ~{measured_rate:.0f} rec/s. "
              f"Extrapolated full-run estimate: ~{est_hours:.2f} hours for {total_count:,} records", flush=True)


if __name__ == "__main__":
    main()
