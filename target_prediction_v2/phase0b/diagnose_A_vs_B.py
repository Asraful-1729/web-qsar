"""
P0b-11: diagnose why Set B (clinical) consistently outperforms Set A
(approved) -- first flagged in Phase 0, never resolved, and sharpened by
the Rev 6 rerun (the entire H4a pooling-vs-singleton effect turned out to
be significant for Set B and not for Set A under correct leakage control).
Rev 5 named three diagnostics explicitly: max-similarity-to-index
distribution, |T|, and document-holdout behaviour. This script computes
all three plus two more that the already-captured data makes cheap
(neighbour-pool density, query's own scaffold-group size), and -- the part
Rev 5 didn't ask for but that actually answers "why" rather than just
"that they differ" -- checks which of these features predicts per-query
recovery success WITHIN each set, so a between-set difference in that
feature's distribution can be read as an actual causal candidate rather
than a coincidental correlation.

Reuses capture_{A,B}.jsonl (already has max_similarity_to_index, the full
neighbour pool, and bestsim_results per query -- no re-querying the index)
plus one direct lookup against the full v1 index for annotated-target
reference depth and query scaffold-group size (fast, in-memory, same
pattern as capture_query_potency.py), plus mechanism_support.json for
document counts (P0b-12's pull, already done).

Usage: python3 diagnose_A_vs_B.py
"""
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PHASE0_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phase0")
sys.path.insert(0, PHASE0_DIR)
sys.path.insert(0, os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "code"))
os.environ.setdefault("TARGET_FISHING_INDEX_DIR", os.path.abspath(
    os.path.join(PHASE0_DIR, "..", "..", "target_fishing_v1_freeze", "target_fishing_index")))

import target_fishing as TF  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")


def load_capture(setname):
    with open(os.path.join(RESULTS_DIR, f"capture_{setname}.jsonl")) as f:
        return [json.loads(line) for line in f]


def scaffold_key(smiles, scaffold):
    if isinstance(scaffold, str) and scaffold:
        return scaffold
    return f"__no_ring__:{smiles}"


def main():
    print("Loading v1 full index for reference-depth and scaffold lookups...", flush=True)
    _, _, full_df = TF._load()
    depth_by_target = full_df.drop_duplicates(subset=["smiles", "target_chembl"]).groupby("target_chembl").size().to_dict()
    scaffold_by_smiles = dict(zip(full_df["smiles"], full_df["murcko_scaffold"]))
    scaffold_groups = full_df.assign(
        scaffold_key=[scaffold_key(s, sc) for s, sc in zip(full_df["smiles"], full_df["murcko_scaffold"])]
    ).drop_duplicates(subset="smiles").groupby("scaffold_key").size().to_dict()

    with open(os.path.join(HERE, "data", "mechanism_support.json")) as f:
        mech_support = json.load(f)

    report = {}
    per_query_rows = {}
    for setname in ["A", "B"]:
        records = load_capture(setname)
        rows = []
        for r in records:
            annotated = r["annotated_targets"]
            t_count = len(annotated)
            depths = [depth_by_target.get(t, 0) for t in annotated]
            mean_depth = statistics.mean(depths) if depths else None
            median_depth = statistics.median(depths) if depths else None
            min_depth = min(depths) if depths else None

            q_scaffold = scaffold_key(r["smiles"], scaffold_by_smiles.get(r["smiles"]))
            own_scaffold_group_size = scaffold_groups.get(q_scaffold, 1)

            neighbours = r["neighbours"]
            n_ge_04 = sum(1 for _, sim, _ in neighbours if sim >= 0.4)
            n_ge_05 = sum(1 for _, sim, _ in neighbours if sim >= 0.5)
            n_ge_07 = sum(1 for _, sim, _ in neighbours if sim >= 0.7)

            mcid = r["molecule_chembl_id"]
            sup = mech_support.get(mcid, {})
            doc_counts = [sup[t]["n_documents"] for t in annotated if t in sup]
            mean_docs = statistics.mean(doc_counts) if doc_counts else None

            # recovery success: is any annotated target in the best_similarity top-10?
            ranked = sorted(r["bestsim_results"], key=lambda ts: -ts[1])
            top10_targets = set(t for t, s in ranked[:10])
            recovered_top10 = bool(set(annotated) & top10_targets)
            recovered_top1 = bool(ranked and ranked[0][0] in annotated)

            rows.append({
                "smiles": r["smiles"],
                "max_similarity_to_index": r["max_similarity_to_index"],
                "n_annotated_targets": t_count,
                "mean_target_reference_depth": mean_depth,
                "median_target_reference_depth": median_depth,
                "min_target_reference_depth": min_depth,
                "own_scaffold_group_size": own_scaffold_group_size,
                "n_neighbours_ge_0.4": n_ge_04,
                "n_neighbours_ge_0.5": n_ge_05,
                "n_neighbours_ge_0.7": n_ge_07,
                "mean_documents_per_annotated_target": mean_docs,
                "recovered_top10": recovered_top10,
                "recovered_top1": recovered_top1,
            })
        per_query_rows[setname] = rows

        def summ(key):
            vals = [row[key] for row in rows if row[key] is not None]
            return {
                "n": len(vals),
                "mean": round(statistics.mean(vals), 3) if vals else None,
                "median": round(statistics.median(vals), 3) if vals else None,
                "p25": round(sorted(vals)[len(vals)//4], 3) if vals else None,
                "p75": round(sorted(vals)[3*len(vals)//4], 3) if vals else None,
            }

        report[setname] = {
            "n_queries": len(rows),
            "max_similarity_to_index": summ("max_similarity_to_index"),
            "n_annotated_targets": summ("n_annotated_targets"),
            "mean_target_reference_depth": summ("mean_target_reference_depth"),
            "median_target_reference_depth": summ("median_target_reference_depth"),
            "min_target_reference_depth": summ("min_target_reference_depth"),
            "own_scaffold_group_size": summ("own_scaffold_group_size"),
            "n_neighbours_ge_0.4": summ("n_neighbours_ge_0.4"),
            "n_neighbours_ge_0.5": summ("n_neighbours_ge_0.5"),
            "n_neighbours_ge_0.7": summ("n_neighbours_ge_0.7"),
            "mean_documents_per_annotated_target": summ("mean_documents_per_annotated_target"),
            "pct_recovered_top10": round(100 * sum(1 for row in rows if row["recovered_top10"]) / len(rows), 1),
            "pct_recovered_top1": round(100 * sum(1 for row in rows if row["recovered_top1"]) / len(rows), 1),
        }

    # within-set correlation: does each feature actually predict recovery, or just
    # differ between sets coincidentally? Compare mean feature value for
    # recovered-top10 vs not-recovered-top10 queries, WITHIN each set.
    predictive = {}
    for setname in ["A", "B"]:
        rows = per_query_rows[setname]
        predictive[setname] = {}
        for key in ["max_similarity_to_index", "n_neighbours_ge_0.4", "n_neighbours_ge_0.5",
                    "mean_target_reference_depth", "min_target_reference_depth", "own_scaffold_group_size"]:
            rec = [row[key] for row in rows if row["recovered_top10"] and row[key] is not None]
            not_rec = [row[key] for row in rows if not row["recovered_top10"] and row[key] is not None]
            predictive[setname][key] = {
                "mean_when_recovered": round(statistics.mean(rec), 3) if rec else None,
                "mean_when_not_recovered": round(statistics.mean(not_rec), 3) if not_rec else None,
                "n_recovered": len(rec), "n_not_recovered": len(not_rec),
            }

    out = {"summary_A_vs_B": report, "within_set_predictive_check": predictive}
    out_path = os.path.join(RESULTS_DIR, "diagnose_A_vs_B_report.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
