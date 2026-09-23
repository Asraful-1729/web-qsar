#!/bin/bash
# Generic self-healing wrapper for any stage script using fetch_common.py's
# checkpointing -- ChEMBL's API has been intermittently returning HTTP 500
# broadly (confirmed not specific to any one query). Keeps relaunching the
# given stage script (resumable, checkpointed) with a cooldown between
# attempts until the checkpoint reaches the known target page count, or a
# max attempt count is hit.
#
# Usage: ./run_stage_until_done.sh <script.py> <checkpoint_file.json> <log_file> <target_pages> [extra script args...]
set -u
cd "$(dirname "$0")"
SCRIPT="$1"
CKPT="$2"
LOG="$3"
TARGET_PAGES="$4"
shift 4
EXTRA_ARGS="$@"
MAX_ATTEMPTS=40

for i in $(seq 1 $MAX_ATTEMPTS); do
    n_done=$(python3 -c "import json
try:
    print(len(json.load(open('$CKPT'))))
except Exception:
    print(0)" 2>/dev/null || echo 0)
    echo "[wrapper attempt $i/$MAX_ATTEMPTS] $CKPT: $n_done/$TARGET_PAGES pages done so far"

    if [ "$n_done" -ge "$TARGET_PAGES" ]; then
        echo "[wrapper] Target reached."
        exit 0
    fi

    python3 "$SCRIPT" $EXTRA_ARGS >> "$LOG" 2>&1

    n_done_after=$(python3 -c "import json
try:
    print(len(json.load(open('$CKPT'))))
except Exception:
    print(0)" 2>/dev/null || echo 0)

    if [ "$n_done_after" -ge "$TARGET_PAGES" ]; then
        echo "[wrapper] Target reached after attempt $i."
        exit 0
    fi

    echo "[wrapper] $n_done_after/$TARGET_PAGES after this attempt, cooling down 45s before retry..."
    sleep 45
done
echo "[wrapper] Reached max attempts ($MAX_ATTEMPTS) without finishing. Manual check needed."
exit 1
