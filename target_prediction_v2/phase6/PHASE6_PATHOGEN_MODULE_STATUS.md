# Phase 6 — Pathogen Module: Real Scoping, Not a Build

**Status: scoped with real ChEMBL/NCBI data, not built.** `BUILD_PLAN.md` explicitly marks this "post-release, off the critical path" and "pocket-embedding models experimental at most" — this pass measures the real landscape so any future work starts from evidence, not assumption, consistent with how every other phase in this program was approached. No production pipeline was built; that would be a genuinely separate, much larger project.

## Bacterial target universe — real numbers, one real bug found and fixed

Built via the same rigor as the human/orthologue target work: fetched the full 11,055-target cross-organism single-protein list, classified each of 700 distinct organism strings via **NCBI's taxonomy API** (not a hardcoded genus list).

**Bug found and fixed**: the first classification pass returned 0 bacterial organisms — NCBI's `esummary` taxonomy response has no `lineage` field (confirmed empirically; that field doesn't exist in this API version), so a lineage-substring check silently failed for everything, including E. coli. Fixed: the real field is `genbankdivision`, directly `"Bacteria"` for bacterial organisms.

**Result: 665 bacterial single-protein targets** (6.0% of the 11,055-target cross-organism universe), across 145 distinct bacterial organism strings.

## Activity data volume — confirms the plan's own "off the critical path" framing, quantified

| | Targets | Base activities (confidence≥8, IC50/Ki/Kd/EC50) | Avg per target |
|---|---:|---:|---:|
| Human tier (Phase 2) | 5,869 | 2,505,626 | ~427 |
| Bacterial tier (this scoping) | 665 | 24,897 | ~37 |

A sample of 12 major WHO-priority pathogens individually (*E. coli, S. aureus, M. tuberculosis, P. aeruginosa, K. pneumoniae, A. baumannii, E. faecium, N. gonorrhoeae, S. pneumoniae, H. pylori, C. difficile, S. typhimurium*) showed 0-86 targets and 0-6,687 activities each — several pathogens (notably *Salmonella typhimurium* in this specific query) return **zero** qualifying targets.

**This directly validates why the plan calls for a mechanism-class fallback**: bacterial query density would almost never reach the density≥8 pooling threshold (`phase3/PHASE3A_DENSITY_REFIT.md`) — the vast majority of bacterial predictions would land in the significantly-negative low-density band (0-1) found during that re-fit. A direct port of the human-tier retrieval pipeline would perform poorly by construction, not by a fixable bug.

## Human-homology proxy — reused the orthologue tier's own method, real result

Same gene-symbol-matching technique already validated for the orthologue tier (`phase2/stage5a_build_orthologue_map.py`): does a bacterial target's gene symbol also appear among the 5,869 human single-protein targets?

**Result: 12/665 bacterial targets (1.8%) share a gene symbol with a human target** — e.g., `inhA` (*M. tuberculosis* enoyl-ACP reductase), `relA`, `fas`. A low match rate, favorable for antibacterial target selection (most bacterial targets are genuinely gene-symbol-distinct from human proteins by this proxy) — though, same caveat the orthologue tier carries: gene-symbol identity is a proxy for homology, not a rigorous sequence-based homology call (e.g., BLAST/OrthoDB), which would need real bioinformatics infrastructure not built here.

## Essentiality annotation — not attempted, genuinely needs external data

Gene essentiality (whether a target is required for pathogen survival — a core druggability signal) is **not derivable from ChEMBL at all**. It requires an external resource (e.g., DEG — Database of Essential Genes, or OGEE) that this session has no programmatic access path to and did not attempt to build. Disclosed as a real gap, not approximated.

## Mechanism-class fallback — a real, unresolved design question, one path identified

ChEMBL's own curated `protein_classifications` field is **sparse-to-absent for bacterial targets** (confirmed: null for `CHEMBL1849`/*M. tuberculosis* InhA, a well-studied first-line TB drug target — if a heavily-studied target has no classification, the tree isn't reliably populated for this tier at all). `EC_NUMBER` synonyms (e.g., "3.2.1.20") **are** present and could serve as an alternative, hierarchical mechanism-class basis (grouping by EC number prefix) — a real, identified path, **not built or tested** in this pass.

## Pocket-embedding models — correctly not attempted

Explicitly "experimental at most" per the plan's own framing — the lowest-priority item in the lowest-priority phase. No work undertaken, consistent with the plan's own prioritization.

## Bottom line

Real numbers now exist for every Phase 6 sub-item except essentiality (external data, not attempted) and pocket-embedding (correctly deprioritized). The bacterial data landscape is confirmed genuinely sparse — any future work here should design around mechanism-class grouping (EC-number-based, not ChEMBL's own classification tree) from the start rather than porting the human-tier pipeline unchanged. This remains correctly off the critical path for v2's release; nothing here should block or delay G1/Phase 7.
