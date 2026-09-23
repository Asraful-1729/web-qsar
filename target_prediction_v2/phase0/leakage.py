"""
Phase 0 (Section 4, point 2) -- per-query leakage control. Removes the query
compound and its near-duplicates (salts/stereoisomers/close analogues) from
the reference BEFORE searching, using TWO independent fingerprints (Morgan,
the model's own; and a path-based topological one, precomputed in
build_leakage_fp2.py) so a duplicate that folds/collides similarly under one
does not slip through silently. Both thresholds the review asked for (0.95
strict, 0.90 looser) are supported; a compound counts as leaked if EITHER
fingerprint reports Tanimoto >= cutoff (union, not intersection -- a stricter
leakage definition, deliberately: this is a benchmark leakage check, not the
production index, so erring toward removing too much is the safe direction).
"""
import numpy as np

_POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint16)


def tanimoto_to_all(q_packed, fps_arr, pop_counts, q_pop=None):
    if q_pop is None:
        q_pop = int(_POPCOUNT_TABLE[q_packed].sum())
    inter = _POPCOUNT_TABLE[np.bitwise_and(fps_arr, q_packed)].sum(axis=1)
    union = q_pop + pop_counts - inter
    return np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)


class LeakageIndex:
    """Wraps the full v1 reference (Morgan fps, aligned df) plus the
       second-fingerprint array (aligned by SMILES, built once by
       build_leakage_fp2.py) so per-query leakage masks can be computed
       against BOTH without re-loading anything per query."""

    def __init__(self, full_fps, full_pop, full_df, fp2_smiles, fp2_arr):
        self.full_fps = full_fps
        self.full_pop = full_pop
        self.full_df = full_df
        # distinct-compound view aligned to fp2's own order
        self.fp2_smiles = fp2_smiles
        self.fp2_arr = fp2_arr
        self.fp2_pop = _POPCOUNT_TABLE[fp2_arr].sum(axis=1)
        self._fp2_index_by_smiles = {s: i for i, s in enumerate(fp2_smiles)}
        # one representative Morgan fingerprint row (into full_fps) per
        # distinct compound, in the SAME order as fp2_smiles -- both are
        # derived from an identical drop_duplicates(subset="smiles") call
        # against the same full_df (see build_leakage_fp2.py), so a plain
        # positional vstack/take, not a per-element Python lookup, aligns
        # them correctly and avoids an 857k-iteration Python loop.
        first_idx = full_df.reset_index().drop_duplicates(subset="smiles", keep="first").set_index("smiles")["index"]
        self._morgan_row_by_smiles = first_idx
        distinct_rows = first_idx.loc[fp2_smiles].to_numpy()
        self._distinct_morgan_fps = full_fps[distinct_rows]
        self._distinct_morgan_pop = full_pop[distinct_rows]

    def query_fingerprints(self, smiles):
        """Returns (morgan_packed, morgan_pop, topo_packed, topo_pop) for a
           compound already present in the v1 index (all Phase 0 queries
           are, by construction -- see build_eval_sets.py)."""
        m_row = int(self._morgan_row_by_smiles[smiles])
        m_packed = self.full_fps[m_row]
        m_pop = int(self.full_pop[m_row])
        t_idx = self._fp2_index_by_smiles.get(smiles)
        if t_idx is None:
            return m_packed, m_pop, None, None
        return m_packed, m_pop, self.fp2_arr[t_idx], int(self.fp2_pop[t_idx])

    def leaked_compound_mask(self, smiles, cutoff):
        """Boolean mask over DISTINCT compounds (fp2's own ordering) that
           count as leaked (query itself, or Tanimoto >= cutoff under
           EITHER fingerprint)."""
        m_packed, m_pop, t_packed, t_pop = self.query_fingerprints(smiles)
        morgan_tanimoto = tanimoto_to_all(m_packed, self._distinct_morgan_fps, self._distinct_morgan_pop, m_pop)
        mask = morgan_tanimoto >= cutoff
        if t_packed is not None:
            topo_tanimoto = tanimoto_to_all(t_packed, self.fp2_arr, self.fp2_pop, t_pop)
            mask = mask | (topo_tanimoto >= cutoff)
        # always remove the exact query compound itself regardless of cutoff
        mask = mask | (self.fp2_smiles == smiles)
        return mask

    def reduced_reference(self, smiles, cutoff):
        """Returns (reduced_fps, reduced_pop, reduced_df) with every row
           belonging to a leaked compound removed."""
        leaked_mask = self.leaked_compound_mask(smiles, cutoff)
        leaked_smiles = set(self.fp2_smiles[leaked_mask].tolist())
        row_mask = ~self.full_df["smiles"].isin(leaked_smiles).values
        return self.full_fps[row_mask], self.full_pop[row_mask], self.full_df.loc[row_mask].reset_index(drop=True)
