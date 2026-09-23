"""
Phase 2 rebuild, Stage 0 (PHASE2_REBUILD_EXECUTION_PLAN.md): pin the exact
ChEMBL release this rebuild runs against (closing the long-standing gap
disclosed in PHASE1_BENCHMARK_CARD.md Sections 1/10 -- this program has
never pinned a release before), and rebuild single_protein_target_ids.json
fresh (v1's copy of this file was a build-time intermediate, not preserved
in the freeze snapshot -- confirmed absent from the repo, so this is a
real rebuild, not a "verify and reuse").

Usage: python3 stage0_pin_release.py
"""
import json
import os
import time
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")


def _get(url, retries=5):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def pin_release():
    status = _get(f"{BASE}/status.json")
    release = {
        "chembl_db_version": status["chembl_db_version"],
        "chembl_release_date": status["chembl_release_date"],
        "pinned_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "status_snapshot": status,
    }
    print(f"Pinned release: {release['chembl_db_version']} ({release['chembl_release_date']})", flush=True)
    return release


def build_single_protein_target_ids():
    ids = []
    url = f"{BASE}/target.json?organism=Homo+sapiens&target_type=SINGLE+PROTEIN&limit=1000"
    t0 = time.time()
    while url:
        d = _get(url)
        ids.extend(t["target_chembl_id"] for t in d["targets"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
        print(f"  {len(ids)} single-protein human targets so far ({time.time()-t0:.0f}s)", flush=True)
    return sorted(set(ids))


def main():
    os.makedirs(DATA_DIR, exist_ok=True)

    release = pin_release()
    with open(os.path.join(DATA_DIR, "release_manifest.json"), "w") as f:
        json.dump(release, f, indent=2)

    print("Building single_protein_target_ids.json fresh (not found preserved from v1's build)...", flush=True)
    ids = build_single_protein_target_ids()
    out_path = os.path.join(DATA_DIR, "single_protein_target_ids.json")
    with open(out_path, "w") as f:
        json.dump(ids, f)
    print(f"Wrote {len(ids)} target ids -> {out_path}", flush=True)
    print(f"(v1's freeze docstring recorded ~5,869 ids as of its own ChEMBL release -- "
          f"a materially different count here would itself be a real, worth-noting finding.)", flush=True)


if __name__ == "__main__":
    main()
