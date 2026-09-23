"""
Phase 6 (Pathogen module) scoping: real bacterial single-protein target
universe from ChEMBL, classified via NCBI taxonomy (not a hardcoded genus
list -- more accurate, and not meaningfully slower given the number of
DISTINCT organism strings across 11,055 single-protein targets is small).

Step 1: re-fetch the full cross-organism single-protein target list
(same pull as phase2/stage5a_build_orthologue_map.py, but this time
saving the raw list -- that script only kept the derived gene-symbol map).
Step 2: extract distinct organism strings, classify each via NCBI
E-utilities (esearch on the taxonomy db + esummary for lineage), checking
for "Bacteria" in the lineage.
Step 3: report real target/activity volume for the resulting bacterial
target set.

Usage: python3 build_bacterial_target_universe.py
"""
import json
import os
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(url, retries=9):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(min(2 ** attempt, 30))


def fetch_all_single_protein_targets():
    out_path = os.path.join(DATA_DIR, "all_single_protein_targets_raw.json")
    if os.path.exists(out_path):
        with open(out_path) as f:
            return json.load(f)
    url = f"{CHEMBL_BASE}/target.json?target_type=SINGLE+PROTEIN&limit=1000"
    targets = []
    t0 = time.time()
    while url:
        d = _get(url)
        targets.extend({"target_chembl_id": t["target_chembl_id"], "organism": t.get("organism"),
                         "pref_name": t.get("pref_name")} for t in d["targets"])
        nxt = d["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
        print(f"  {len(targets)} targets fetched ({time.time()-t0:.0f}s)", flush=True)
    with open(out_path, "w") as f:
        json.dump(targets, f)
    return targets


def classify_organism_ncbi(organism_name, retries=5):
    """Returns True if the organism's NCBI genbankdivision is 'Bacteria',
    False if resolved but not bacterial, None if unresolvable.

    FOUND AND FIXED, not assumed: esummary's taxonomy response has no
    'lineage' field (confirmed live) -- that field simply doesn't exist in
    this NCBI API version, so checking for it always returned "" and every
    organism silently classified as non-bacterial (0/700 bacterial on the
    first run, including E. coli). The real field is 'genbankdivision',
    directly set to 'Bacteria' for bacterial organisms -- cleaner than a
    lineage substring match anyway."""
    try:
        q = urllib.parse.quote(f'"{organism_name}"[Scientific Name]')
        url = f"{EUTILS}/esearch.fcgi?db=taxonomy&term={q}&retmode=json"
        d = _get(url, retries=retries)
        ids = d.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None
        taxid = ids[0]
        url2 = f"{EUTILS}/esummary.fcgi?db=taxonomy&id={taxid}&retmode=json"
        d2 = _get(url2, retries=retries)
        division = d2.get("result", {}).get(taxid, {}).get("genbankdivision", "")
        return division == "Bacteria"
    except Exception:
        return None


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    print("Fetching full cross-organism single-protein target list...", flush=True)
    targets = fetch_all_single_protein_targets()
    print(f"{len(targets)} single-protein targets total", flush=True)

    organisms = sorted(set(t["organism"] for t in targets if t.get("organism")))
    print(f"{len(organisms)} distinct organism strings to classify", flush=True)

    class_path = os.path.join(DATA_DIR, "organism_bacteria_classification.json")
    classification = {}
    if os.path.exists(class_path):
        with open(class_path) as f:
            classification = json.load(f)
        print(f"Resuming: {len(classification)} already classified", flush=True)

    remaining = [o for o in organisms if o not in classification]
    t0 = time.time()
    for i, org in enumerate(remaining):
        classification[org] = classify_organism_ncbi(org)
        if (i + 1) % 50 == 0:
            with open(class_path, "w") as f:
                json.dump(classification, f)
            print(f"  {i+1}/{len(remaining)} classified ({time.time()-t0:.0f}s)", flush=True)
    with open(class_path, "w") as f:
        json.dump(classification, f)

    n_bacterial = sum(1 for v in classification.values() if v is True)
    n_unresolved = sum(1 for v in classification.values() if v is None)
    print(f"\n{n_bacterial} organisms classified as Bacteria, "
          f"{n_unresolved} unresolvable, "
          f"{len(classification) - n_bacterial - n_unresolved} classified non-bacterial", flush=True)

    bacterial_targets = [t for t in targets if classification.get(t.get("organism")) is True]
    print(f"\n{len(bacterial_targets)} bacterial single-protein targets found", flush=True)

    out_path = os.path.join(DATA_DIR, "bacterial_single_protein_targets.json")
    with open(out_path, "w") as f:
        json.dump(bacterial_targets, f, indent=2)
    print(f"Wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
