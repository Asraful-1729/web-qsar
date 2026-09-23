"""
Phase 2 follow-up: document-level split, the last of the "four split
types" G1 asks for. Closer to buildable than temporal turned out to be
(document_chembl_id was already present per raw row, per pair via Stage
7's now-exposed `document_chembl_ids` field -- no new ChEMBL pull needed,
unlike the date backfill).

Partition unit: document_chembl_id groups (a compound is grouped with
every OTHER compound that shares at least one contributing document --
transitively, via union-find, since documents and compounds form a
bipartite graph and a naive per-document split would let a multi-document
compound straddle train/test). This is the real leakage concern rev5
flagged: an entire assay campaign (one document, many compounds/targets)
should not have some rows in train and others in test.

FOUND AND FIXED, not assumed: a naive union-find over ALL documents
produces one catastrophic cluster (1,016,961 of 1,347,039 compounds, 75%)
because Stage 3's mega-document (CHEMBL1201862, the same one found to
dominate the temporal-holdout gap) alone connects nearly a million
compounds -- a bulk-deposit artifact, not a genuine shared-assay-campaign
leakage risk. First run of this script produced an 86%/14% split instead
of 80/20 because of it. FIX: documents connecting more than
MEGA_DOCUMENT_THRESHOLD compounds are excluded from the grouping logic
entirely (their compounds still get split normally via their OTHER,
smaller-scale document connections) -- documented here, not silently
patched.

Usage: python3 build_document_split.py
"""
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
TEST_FRACTION = 0.20
SEED = 42
MEGA_DOCUMENT_THRESHOLD = 1000  # documents connecting more compounds than this are excluded from grouping


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def main():
    print("Pass 1: loading pairs, counting per-document compound degree...", flush=True)
    compound_docs = {}  # smiles -> set of doc ids (for reporting)
    doc_compounds = {}  # doc id -> set of compound smiles (to find mega-documents)
    n_pairs = 0
    n_pairs_with_doc = 0
    with open(os.path.join(DATA_DIR, "stage7_native_human_pairs.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            n_pairs += 1
            smi = r["smiles"]
            docs = r.get("document_chembl_ids") or []
            if not docs:
                continue
            n_pairs_with_doc += 1
            compound_docs.setdefault(smi, set()).update(docs)
            for d in docs:
                doc_compounds.setdefault(d, set()).add(smi)

    print(f"{n_pairs} pairs total, {n_pairs_with_doc} with >=1 document ({100*n_pairs_with_doc/n_pairs:.1f}%)", flush=True)

    mega_docs = {d for d, cs in doc_compounds.items() if len(cs) > MEGA_DOCUMENT_THRESHOLD}
    print(f"{len(mega_docs)} mega-document(s) excluded from grouping "
          f"(connect >{MEGA_DOCUMENT_THRESHOLD} compounds each): "
          f"{sorted(mega_docs, key=lambda d: -len(doc_compounds[d]))[:5]} "
          f"(largest connects {max((len(doc_compounds[d]) for d in mega_docs), default=0)} compounds)", flush=True)

    print("Pass 2: building union-find over non-mega documents...", flush=True)
    uf = UnionFind()
    for smi, docs in compound_docs.items():
        cnode = ("C", smi)
        uf.find(cnode)
        for d in docs:
            if d in mega_docs:
                continue
            uf.union(cnode, ("D", d))

    # group compounds by their union-find root (their connected document-sharing cluster)
    groups = {}
    for smi in compound_docs:
        root = uf.find(("C", smi))
        groups.setdefault(root, []).append(smi)

    sizes = sorted((len(v) for v in groups.values()), reverse=True)
    print(f"{len(groups)} distinct document-connected compound clusters "
          f"(largest: {sizes[0]}, median: {sizes[len(sizes)//2]}, "
          f"n singletons: {sum(1 for s in sizes if s == 1)})", flush=True)

    total = sum(len(v) for v in groups.values())
    target_test_n = int(round(total * TEST_FRACTION))

    # size-aware assignment, not pure random shuffle: even after excluding
    # direct mega-documents, transitive chaining through moderate-degree
    # documents still produced one large residual cluster (580,958
    # compounds, confirmed by the first run of this fix) that a random
    # shuffle would let dominate whichever side it landed on, producing a
    # ~50/50 split instead of 80/20. Any cluster too large to fit within
    # the test budget is forced to train (it can never fit a 20% test set
    # regardless of shuffle order) -- remaining clusters are randomly
    # assigned to fill the target ratio.
    import random
    rng = random.Random(SEED)
    oversized = [gk for gk in groups if len(groups[gk]) > target_test_n]
    fittable = [gk for gk in groups if len(groups[gk]) <= target_test_n]
    rng.shuffle(fittable)
    print(f"{len(oversized)} cluster(s) too large for the {TEST_FRACTION:.0%} test budget "
          f"(target={target_test_n}) -- forced to train", flush=True)

    group_keys = fittable  # oversized go straight to train below
    test_groups, train_groups = [], list(oversized)
    running = 0
    for gk in group_keys:
        if running < target_test_n:
            test_groups.append(gk)
            running += len(groups[gk])
        else:
            train_groups.append(gk)

    test_smiles = [s for gk in test_groups for s in groups[gk]]
    train_smiles = [s for gk in train_groups for s in groups[gk]]

    out = {
        "method": "document-connected compound clusters (union-find over shared document_chembl_id), "
                  "random assignment of whole clusters to test/train -- no cluster straddles the split",
        "seed": SEED,
        "n_pairs_total": n_pairs,
        "n_pairs_with_document": n_pairs_with_doc,
        "pct_pairs_with_document": round(100 * n_pairs_with_doc / n_pairs, 1),
        "n_compounds_with_document": len(compound_docs),
        "n_document_connected_clusters": len(groups),
        "largest_cluster_size": sizes[0],
        "document_test": {"n_compounds": len(test_smiles), "n_clusters": len(test_groups), "smiles": test_smiles},
        "document_train": {"n_compounds": len(train_smiles), "n_clusters": len(train_groups)},
    }
    out_path = os.path.join(DATA_DIR, "document_split.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nDocument test: {len(test_smiles)} compounds, {len(test_groups)} clusters", flush=True)
    print(f"Document train: {len(train_smiles)} compounds, {len(train_groups)} clusters", flush=True)
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
