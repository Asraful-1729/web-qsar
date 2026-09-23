"""
Phase 2 rebuild -- shared concurrent-paginated-fetch-with-checkpoint core,
extracted from Stage 1's script (stage1_fetch_base_evidence.py) after that
script's design was validated at full production scale (885 rec/s, 30
concurrent workers, 0 data-integrity issues on a 100K-record pilot).
Stage 1's own script is left as-is (already running/validated) rather than
refactored onto this module mid-flight; Stages 2-4 use this shared core
directly instead of re-duplicating it.

Every stage using this module gets, for free: 30-worker concurrent paged
fetching, per-page-offset checkpointing (safe under out-of-order
concurrent completion), resume-from-partial, and a --pilot cap for a
measured-throughput test run before committing to the full pull -- the
same properties Stage 1 validated, not a new untested pattern.
"""
import json
import math
import os
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

PAGE_LIMIT = 1000
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


def run_concurrent_fetch(base_query_url, out_path, ckpt_path, pilot=None,
                          workers=N_WORKERS, row_transform=None):
    """base_query_url: full activity.json query string, WITHOUT limit/offset.
    row_transform(record) -> record: applied to every row before writing
        (e.g. add is_censored=True, compute pchembl_value locally). May
        return None to drop a row (e.g. a unit/value that fails validation).
    Returns (n_pages_done_this_run, elapsed_seconds, total_count)."""

    def fetch_page(offset):
        url = f"{base_query_url}&limit={PAGE_LIMIT}&offset={offset}"
        d = _get(url)
        return offset, d["activities"]

    probe = _get(f"{base_query_url}&limit=1")
    total_count = probe["page_meta"]["total_count"]
    print(f"Total matching records available: {total_count}", flush=True)

    n_pages = math.ceil(total_count / PAGE_LIMIT)
    all_offsets = [i * PAGE_LIMIT for i in range(n_pages)]
    if pilot:
        n_pilot_pages = math.ceil(pilot / PAGE_LIMIT)
        all_offsets = all_offsets[:n_pilot_pages]
        print(f"PILOT MODE: capping at {len(all_offsets)} pages (~{pilot} records)", flush=True)

    completed = set()
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            completed = set(json.load(f))
    remaining = [o for o in all_offsets if o not in completed]
    if completed:
        print(f"Resuming: {len(completed)}/{len(all_offsets)} pages already done, "
              f"{len(remaining)} remaining", flush=True)

    def save_ckpt():
        with open(ckpt_path, "w") as f:
            json.dump(sorted(completed), f)

    t0 = time.time()
    write_lock = threading.Lock()
    pages_done_this_run = 0
    n_dropped = 0

    with open(out_path, "a") as f:
        with ThreadPoolExecutor(max_workers=workers) as ex:
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
                        if row_transform:
                            r2 = row_transform(r)
                            if r2 is None:
                                n_dropped += 1
                                continue
                            r = r2
                        f.write(json.dumps(r) + "\n")
                    f.flush()
                    completed.add(off)
                    pages_done_this_run += 1
                    if pages_done_this_run % 10 == 0 or len(completed) == len(all_offsets):
                        save_ckpt()

                elapsed = time.time() - t0
                rec_rate = pages_done_this_run * PAGE_LIMIT / elapsed if elapsed > 0 else 0
                print(f"  {len(completed)}/{len(all_offsets)} pages "
                      f"({elapsed:.0f}s elapsed this run, ~{rec_rate:.0f} rec/s, "
                      f"{n_dropped} dropped by row_transform)", flush=True)

    save_ckpt()
    elapsed = time.time() - t0
    print(f"\nDone: {len(completed)}/{len(all_offsets)} pages -> {out_path} "
          f"({elapsed:.0f}s, {n_dropped} rows dropped)", flush=True)
    if pilot and elapsed > 0 and pages_done_this_run > 0:
        measured_rate = pages_done_this_run * PAGE_LIMIT / elapsed
        est_min = total_count / measured_rate / 60
        print(f"Measured rate this run: ~{measured_rate:.0f} rec/s. "
              f"Extrapolated full-run estimate: ~{est_min:.1f} min for {total_count:,} records", flush=True)

    return pages_done_this_run, elapsed, total_count
