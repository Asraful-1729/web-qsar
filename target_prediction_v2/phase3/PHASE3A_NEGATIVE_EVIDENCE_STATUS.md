# Phase 3A — Negative-Evidence Term: Infrastructure Gap, Not Yet Testable

**Status: BLOCKED on a schema extension, disclosed rather than approximated.**

The current `capture_scaffold_strict_v2.py` neighbour schema stores `[target_chembl, pchembl_or_null]` per neighbour-target link. A `pchembl=null` entry is now genuinely ambiguous between three different real cases: no evidence at all, evidence with no computed pchembl, and — since `v2_index/compounds.csv.gz` includes `is_measured_inactive=True`-only pairs from Stage 7 without carrying that flag through — **explicitly measured-inactive evidence** (item 6's ~490K records). These three cases need different treatment in a negative-evidence scoring term (only the third should count as a penalty signal); collapsing them all to `null` makes a correct test impossible with the current capture data.

**Not approximated or faked.** Fixing this requires: (1) adding `is_measured_inactive` to `v2_index/compounds.csv.gz`'s schema (a small change to `build_v2_index_files.py`), (2) extending `capture_scaffold_strict_v2.py`'s neighbour-tuple format to carry it through, (3) re-running the n=800/set captures a third time, (4) extending `score.py` with an actual negative-evidence scoring term (not yet written — no existing function to reuse here, unlike popularity/potency which already had implementations to test).

**Recommendation**: schedule as its own follow-up pass, same rigor as the other two levers (implement, test at real scale with CI, don't assume it helps). Left open here rather than blocking the rest of Phase 3A/3B-7.
