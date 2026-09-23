"""
Phase 0 (target_prediction_v2_procedure_rev3.md, Section 4) — external reference
data pull. This is NEW data, independent of target_fishing_v1_freeze's own index:
it supplies drug-realistic evaluation SETS A and B and the primary/mechanism-
target ground truth needed for Section 3.1's three-level target definitions.
v1's own index (compounds.csv.gz) only stores standardized SMILES, no
molecule_chembl_id, so it cannot answer "which of these are approved drugs" or
"what is compound X's PRIMARY target" on its own -- both questions need this
separate ChEMBL pull.

Two ChEMBL REST endpoints, both small and cheap (checked interactively before
writing this: 4,225 + 1,892 + 9,054 + 1,613 molecules across max_phase 1-4;
7,561 mechanism records total -- a few thousand paginated requests at most,
not a repeat of the original activities.jsonl-scale pull):

  1. molecule.json?max_phase={1,2,3,4} -> molecule_chembl_id, canonical_smiles,
     max_phase. This is Sections 4.1's Set A (max_phase=4, approved) and the
     basis for Set B (max_phase in {1,2,3}, clinical/preclinical -- see
     build_eval_sets.py's docstring for why "preclinical" proper is NOT
     represented here and how that gap is disclosed rather than papered over).

  2. mechanism.json -> molecule_chembl_id, target_chembl_id, action_type,
     mechanism_of_action, direct_interaction, molecular_mechanism. This is
     ChEMBL's own curated drug-mechanism table -- the "Primary/intended
     target" operationalisation Section 3.1 asks to define, verified here to
     actually exist at usable scale (answers "Decisions needed" item 2 in the
     procedure doc empirically: 7,561 records is enough to build primary-
     target metrics on, though coverage per-drug is checked in
     build_eval_sets.py, not assumed here).

Output (this directory's data/ subfolder):
  - drugs_raw.csv       molecule_chembl_id, canonical_smiles, max_phase
  - mechanisms_raw.csv  molecule_chembl_id, target_chembl_id, action_type,
                        mechanism_of_action, direct_interaction,
                        molecular_mechanism

Usage: python3 fetch_chembl_reference.py
"""
import csv
import json
import os
import time
import urllib.request

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
BASE = "https://www.ebi.ac.uk/chembl/api/data"
PAGE = 1000


def _get(url, retries=5):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return json.load(resp)
        except Exception as e:
            if attempt == retries - 1:
                raise
            wait = 2 ** attempt
            print(f"  retry {attempt+1}/{retries} after {e} (sleeping {wait}s)", flush=True)
            time.sleep(wait)


def fetch_molecules_by_phase(max_phase):
    rows = []
    url = (f"{BASE}/molecule.json?max_phase={max_phase}&limit={PAGE}"
           f"&only=molecule_chembl_id,max_phase,molecule_structures")
    t0 = time.time()
    n = 0
    while url:
        d = _get(url)
        for m in d["molecules"]:
            struct = m.get("molecule_structures") or {}
            smi = struct.get("canonical_smiles")
            if not smi:
                continue
            rows.append({
                "molecule_chembl_id": m["molecule_chembl_id"],
                "canonical_smiles": smi,
                "max_phase": max_phase,
            })
        n += len(d["molecules"])
        nxt = d["page_meta"].get("next")
        url = (f"https://www.ebi.ac.uk{nxt}") if nxt else None
        if n % 2000 == 0 or url is None:
            print(f"  max_phase={max_phase}: {n}/{d['page_meta']['total_count']} ({time.time()-t0:.0f}s)", flush=True)
    return rows


def fetch_mechanisms():
    rows = []
    url = (f"{BASE}/mechanism.json?limit={PAGE}"
           f"&only=molecule_chembl_id,target_chembl_id,action_type,mechanism_of_action,"
           f"direct_interaction,molecular_mechanism")
    t0 = time.time()
    n = 0
    while url:
        d = _get(url)
        for m in d["mechanisms"]:
            if not m.get("molecule_chembl_id") or not m.get("target_chembl_id"):
                continue
            rows.append({
                "molecule_chembl_id": m["molecule_chembl_id"],
                "target_chembl_id": m["target_chembl_id"],
                "action_type": m.get("action_type"),
                "mechanism_of_action": m.get("mechanism_of_action"),
                "direct_interaction": m.get("direct_interaction"),
                "molecular_mechanism": m.get("molecular_mechanism"),
            })
        n += len(d["mechanisms"])
        nxt = d["page_meta"].get("next")
        url = (f"https://www.ebi.ac.uk{nxt}") if nxt else None
        if n % 2000 == 0 or url is None:
            print(f"  mechanisms: {n}/{d['page_meta']['total_count']} ({time.time()-t0:.0f}s)", flush=True)
    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Fetching molecules by max_phase (1, 2, 3, 4)...", flush=True)
    all_drugs = []
    for mp in (4, 3, 2, 1):
        all_drugs.extend(fetch_molecules_by_phase(mp))
    drugs_path = os.path.join(OUT_DIR, "drugs_raw.csv")
    with open(drugs_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["molecule_chembl_id", "canonical_smiles", "max_phase"])
        w.writeheader()
        w.writerows(all_drugs)
    print(f"  wrote {len(all_drugs)} rows -> {drugs_path}", flush=True)

    print("Fetching mechanism.json (drug-mechanism / primary-target table)...", flush=True)
    mechs = fetch_mechanisms()
    mech_path = os.path.join(OUT_DIR, "mechanisms_raw.csv")
    with open(mech_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["molecule_chembl_id", "target_chembl_id", "action_type",
                                           "mechanism_of_action", "direct_interaction", "molecular_mechanism"])
        w.writeheader()
        w.writerows(mechs)
    print(f"  wrote {len(mechs)} rows -> {mech_path}", flush=True)


if __name__ == "__main__":
    main()
