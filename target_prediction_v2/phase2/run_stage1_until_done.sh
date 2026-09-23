#!/bin/bash
# Self-healing wrapper: ChEMBL's API has been intermittently returning
# HTTP 500 (confirmed broadly, not specific to our query -- even bare
# status.json/activity.json?limit=1 failed during testing). Keeps
# relaunching stage1_fetch_base_evidence.py (resumable, checkpointed)
# with a cooldown between attempts until all 2506 pages are done or a
# max attempt count is hit.
cd "$(dirname "$0")"
MAX_ATTEMPTS=30
TARGET_PAGES=2506

for i in $(seq 1 $MAX_ATTEMPTS); do
    n_done=$(python3 -c "import json; print(len(json.load(open('data/stage1_completed_offsets.json'))))" 2>/dev/null || echo 0)
    echo "[wrapper attempt $i/$MAX_ATTEMPTS] $n_done/$TARGET_PAGES pages done so far"
    if [ "$n_done" -ge "$TARGET_PAGES" ]; then
        echo "[wrapper] All pages done."
        exit 0
    fi
    python3 stage1_fetch_base_evidence.py >> stage1_full_run.log 2>&1
    n_done=$(python3 -c "import json; print(len(json.load(open('data/stage1_completed_offsets.json'))))" 2>/dev/null || echo 0)
    if [ "$n_done" -ge "$TARGET_PAGES" ]; then
        echo "[wrapper] All pages done after attempt $i."
        exit 0
    fi
    echo "[wrapper] Not done yet ($n_done/$TARGET_PAGES), cooling down 45s before retry..."
    sleep 45
done
echo "[wrapper] Reached max attempts ($MAX_ATTEMPTS) without finishing. Manual check needed."
exit 1
