# Similarity Search — Complete Technical & Scientific Documentation

**Scope:** the natural-product similarity search feature — the Similarity Search tab, its one-time shared index download, and the search engine underneath (`backend/similarity.py`, `backend/scripts/build_similarity_index.py`, the `/api/similarity/*` routes in `backend/app.py`, `frontend/src/tabs/SimilarityTab.tsx`). Written directly from the current source code, not from memory.

**Audience:** anyone who needs to know exactly what "similar" means here, where the reference compound database comes from, and what the returned numbers do and don't imply.

---

## 1. Objective

Given a **known drug or lead compound** (any SMILES — not necessarily a natural product itself), find **structurally similar natural products** from a large reference database, to support "does nature already have something like this" reasoning — a common early step in natural-product-based drug discovery (e.g., starting from a known pharmacophore and asking which real, already-characterized natural products resemble it). This is a pure cheminformatics similarity search — no biological activity, target, or potency claim is made or implied by a result appearing here; that is what the rest of the app's other tabs (Predict, Docking, Target Fishing) are for.

---

## 2. Scientific and methodological basis

### 2.1 The reference database — COCONUT

The search index is built from **[COCONUT](https://coconut.naturalproducts.net/)** (COlleCtion of Open Natural prodUcTs), a real, published, actively-maintained open natural-product database. Specifically, the COCONUT "Lite" CSV export — **CC0-licensed** (public domain — free to use, modify, and redistribute, no attribution legally required, though this app credits it throughout). As of the September 2026 release used to build the current index: **~738,000 compounds**, of which the build successfully parsed and fingerprinted **the number recorded in the index's own `manifest.json`** at build time (a small fraction of raw rows fail to parse with RDKit and are dropped, never silently substituted).

This is the same category of asset as this app's other large reference resources (the target-fishing ChEMBL index, the QSAR model buckets): **built once, offline, and shipped/downloaded as a static artifact** — never computed live from a live database connection.

### 2.2 The similarity method

**Morgan/ECFP4 fingerprints** (radius 2, 2048 bits) — the same fingerprinting convention used consistently everywhere else in this app (QSAR featurization, target fishing, DUD-E-style decoy generation), so a "similarity" number here means the same structural-similarity concept a user would already be familiar with from those other tabs.

**Tanimoto coefficient** on the packed fingerprint bits is the sole ranking criterion — results are sorted purely by descending Tanimoto similarity to the query, nothing else factors into rank.

**Murcko scaffold match** — computed for the query and every candidate independently, and reported as a boolean per result: does this candidate's Bemis–Murcko scaffold (the core ring system, side chains stripped) match the query's own scaffold *exactly* (string-identical canonical SMILES)? This is a real, additional structural signal a user can use to distinguish "similar overall shape/substituents, different core" from "genuinely the same molecular scaffold" — useful context when ranking by Tanimoto alone can't make this distinction.

**Maximum common substructure (MCS)** — computed for the **top `mcs_top_n` results only** (default 10), not the full result set. `rdFMCS.FindMCS` is comparatively expensive per compound pair (it searches for the largest shared substructure, capped by a 2-second per-pair timeout in this implementation) — running it against every one of up to 200 returned hits would be needlessly slow when only the closest matches are typically of real interest. Each of those top results is annotated with the actual shared substructure (`mcs_smarts`) and its atom count (`mcs_n_atoms`) — real, computed values, not estimated from the Tanimoto score.

### 2.3 Performance engineering — why this is fast

The module's own docstring records a concrete, measured before/after: an earlier approach that constructed one RDKit `ExplicitBitVect` object per indexed compound and compared them one at a time measured **~95 seconds just to load the index**. The current implementation instead:

1. Stores every compound's fingerprint as **packed bits** (256 bytes per compound, `uint8` numpy array) rather than as RDKit bit-vector objects.
2. Computes Tanimoto similarity via **bitwise AND + a precomputed popcount lookup table** (`_POPCOUNT_TABLE`, mapping each possible byte value 0–255 to its bit count) applied across the *entire* index array at once with numpy — one vectorized operation against all ~738K compounds simultaneously, not a per-compound Python loop.
3. Precomputes each indexed compound's own total bit-population count once, at index-load time, and reuses it for every subsequent query (rather than recomputing it per search).

Measured result: **well under 1 second to both load the index and run the first query** — and the two implementations were **cross-validated to produce identical Tanimoto scores to 3 decimal places** before the faster one was trusted, i.e. the speed came with a verified-correctness check, not just a reasonable assumption that a bitwise reimplementation would match.

### 2.4 What a result does and doesn't mean

A high Tanimoto score means: this natural product shares a substantial fraction of the same Morgan/ECFP4 circular substructural environments as the query. It does **not** mean the natural product has been shown (or is predicted by this feature) to have the same biological activity — this feature makes no potency, target, or binding claim whatsoever. A researcher would typically follow up a promising similarity hit with this app's other tools (QSAR prediction, docking, target fishing) against a specific biological question, using the similarity hit only as a starting point for "here's a real, already-known natural product worth investigating further."

---

## 3. Architecture — module map

### Backend

| Module | Responsibility |
|---|---|
| `similarity.py` | Everything the live app does: index availability check, the one-time shared-resource download (with progress tracking and cancellation), the packed-bit Tanimoto search engine (§2.2–2.3), Murcko scaffold matching, top-hit MCS computation. |
| `scripts/build_similarity_index.py` | The **offline, one-off** index builder: raw COCONUT CSV → per-compound Morgan fingerprint + Murcko framework + metadata → three output files (`fingerprints.npy`, `compounds.csv`, `manifest.json`). Not run by the live app; run once by a maintainer, and the resulting files are what gets uploaded to remote storage for on-demand download. |
| `app.py` (`/api/similarity/*`) | The public API surface: index status, download start/progress/cancel, and the search endpoint itself. |

### Frontend

| File | Responsibility |
|---|---|
| `tabs/SimilarityTab.tsx` | The tab itself: index-readiness state machine (loading → unavailable/downloading → ready), the one-time download UI, the search form, and an inline results preview. |

---

## 4. Data provenance and index building (offline, maintainer-side)

Not part of the interactive user journey, but part of the same scientific pipeline, documented for completeness and auditability:

1. A raw COCONUT "Lite" CSV export is obtained directly from `coconut.naturalproducts.net/download`.
2. `build_similarity_index.py` streams through it in chunks (default 20,000 rows at a time — the full file is too large to hold naively), and for every row with a parseable `canonical_smiles`:
   - Computes its Morgan/ECFP4 fingerprint (radius 2, 2048 bits), packed to bytes.
   - Re-canonicalizes the SMILES through RDKit (so the stored SMILES is this app's own canonical form, not necessarily COCONUT's).
   - Retains: id, canonical SMILES, name, molecular weight, Murcko framework, NP-classifier pathway/superclass/class, and chemical superclass — everything the search UI needs, without ever needing to re-parse the full ~550MB raw file at query time.
3. Writes three **row-aligned** outputs (`fingerprints.npy` row *i* corresponds to `compounds.csv` row *i* — the search code relies on this alignment directly, with no join key needed at query time): the packed fingerprint array, the compound metadata table, and a manifest recording real counts (rows read vs. successfully parsed), the fingerprint convention used, the source, and a build timestamp.
4. The resulting three files are packaged and uploaded to the same public, credential-free remote storage bucket the per-target model/docking downloads use (`downloads.py`'s `DOWNLOAD_BASE_URL`), under `similarity_index/similarity_index.zip` — from which the live app's own download flow (§5.2) fetches it.

---

## 5. The user journey

### 5.1 Readiness

On opening the tab, `/api/similarity/status` is checked once. Three states follow:
- **Unavailable** — the index hasn't been downloaded to this machine yet. A notice states its real, fixed size (**~87 MB, downloaded once**) and offers a **"Download similarity index"** button.
- **Downloading** — a progress bar (real bytes downloaded / real total bytes from the response's `Content-Length`, in MB) with a **Stop** button.
- **Ready** — the search form appears.

### 5.2 The one-time shared download

This is deliberately **not** routed through the app's general per-target download-gate machinery (`downloads.py`, used for QSAR model buckets and docking data) — the module docstring states why directly: *"This is a SINGLE shared resource (not per-target), fetched on demand via `/api/similarity/download` rather than the per-target download-gate machinery."* One index serves every search, regardless of which target or project a user is working on, so it has its own small, simpler download-job system:

1. A background thread streams the remote ZIP to a temporary file (`.part` suffix) inside the index directory itself, checking a cooperative cancellation flag once per 1 MB chunk — the same cooperative-cancellation pattern used by every other long-running download/job in this app.
2. Once fully downloaded, it's extracted to a temporary subdirectory, and each of the three expected files is moved to a `.new`-suffixed path first, then **atomically renamed** (`os.replace`) over the real file — so a reader that happens to check `available()` mid-extraction never sees a partially-replaced index.
3. The in-memory index cache (`functools.lru_cache` on `_load()`) is cleared once the swap completes, so the very next search picks up the freshly-downloaded data without requiring a server restart.

### 5.3 Search

**Query SMILES** — any valid SMILES string; there is no requirement that the query itself be a natural product (the intended use is a *known drug or lead compound* as the query, per the sidebar's own description).

**Minimum Tanimoto similarity** — a numeric threshold (default 0.4, the same conventional default used elsewhere in this app's similarity-based tools), adjustable in steps of 0.05.

**"Search"** calls `/api/similarity/search`, which runs the full pipeline (§2.2) against the currently-loaded index (up to 200 results, `top_n`, capped server-side) and returns them already sorted by descending Tanimoto.

**Results** are shown as a compact inline list (in the sidebar itself, directly below the search form, rather than in the tab's main content area) — the component's own code comment states this plainly as a deliberate scope decision for a single-query tool, not an oversight: *"the sidebar only shows a compact preview... for a single-query tool a simple top-N preview here is enough."* Each result shows: name (or id, if unnamed), the Tanimoto score as a badge, the full SMILES (truncated, full string on hover), whether it shares the query's exact Murcko scaffold, its chemical superclass (when known), and — for the top `mcs_top_n` results — the MCS atom count.

A summary line above the results states the real counts: how many matches cleared the threshold, out of how many total compounds were searched (e.g. "*12 match(es) above 0.4 similarity, out of 738,204 indexed compounds*").

---

## 6. Background automatic work and decisions — consolidated

### 6.1 Always automatic

- Index-availability check on tab load.
- Fingerprinting the query molecule and computing Tanimoto against the entire index, on every search.
- Murcko scaffold computation and exact-match comparison, for the query and every result.
- MCS computation, automatically limited to the top `mcs_top_n` (10) results — never run against the full result set, invisibly bounding cost regardless of how large `top_n` is set.
- Atomic, safe-to-interrupt index file replacement during a download (§5.2) — a crashed or cancelled download can never leave a half-written index in place that a subsequent search would silently load.
- In-memory index cache invalidation the instant a fresh download completes.

### 6.2 Automatic by default, adjustable in the UI

- Minimum Tanimoto threshold (default 0.4).
- Number of results returned (`top_n`, capped at 200 server-side; not exposed as a separate UI control beyond the implicit default).

### 6.3 Never automatic

- Downloading the index in the first place — always an explicit user click, even though the tab clearly signals when it's needed.
- Submitting a search — always an explicit action on a specific query molecule.

---

## 7. Limitations — stated directly

- **This is a pure structural-similarity tool. It makes no biological, activity, or safety claim about any result.** A high Tanimoto score to a known active drug does not imply the matched natural product shares its activity, potency, or safety profile — only that it shares substantial 2D structural substructure content by the Morgan/ECFP4 convention.
- **2D circular-fingerprint (Morgan/ECFP4) similarity has well-known blind spots** shared with every other tool in this app that uses the same convention: it can miss "activity cliffs" (near-identical fingerprints, very different activity) and can under- or over-weight similarity depending on scaffold size and substituent pattern — a documented, general limitation of the method family, not specific to this implementation.
- **The Murcko scaffold match is an exact string match on the canonical scaffold SMILES.** A candidate with a scaffold that is chemically extremely close but not string-identical (e.g. differing by one ring-fusion atom) will show `scaffold_match: false` — this is a strict, not fuzzy, comparison.
- **MCS is capped by a 2-second per-pair timeout and only computed for the top 10 hits by default** — a genuinely large maximum common substructure search that would need longer than 2 seconds is truncated by RDKit's own timeout behavior (a partial/approximate MCS may be returned rather than the true maximum), and results ranked 11th or lower never get an MCS annotation at all regardless of how interesting they might be.
- **The reference database is a snapshot of COCONUT at whatever release the index was last built from** — new natural products added to COCONUT after that build date will not appear until the index is rebuilt and redownloaded; the manifest's `built_at` timestamp is the authoritative record of the index's actual age, and the UI does not currently surface this timestamp to the user.
- **A compound's absence from a search result does not mean no similar natural product exists** — it means no compound in *this specific COCONUT snapshot* cleared the chosen similarity threshold; COCONUT itself, while large, is not an exhaustive record of all known or all possible natural products.
- **Results render only as a compact inline sidebar preview, not a full sortable/filterable table in the main content area** — a known, deliberate scope simplification for what is currently a single-query tool, recorded directly in the component's own code comments rather than treated as a hidden gap.

---

## 8. API reference (similarity-search-specific endpoints)

| Endpoint | Purpose |
|---|---|
| `GET /api/similarity/status` | Whether the index is downloaded and ready. |
| `POST /api/similarity/download` | Start the one-time index download (or report it's already installed). |
| `GET /api/similarity/download/progress/{job_id}` | Poll download progress (bytes done/total, state). |
| `POST /api/similarity/download/cancel/{job_id}` | Cancel an in-progress download. |
| `POST /api/similarity/search` | Run a similarity search: `smiles`, `threshold` (0–1, default 0.4), `top_n` (1–200, default 50). |

---

## 9. Document provenance

Written by reading, in full: `backend/similarity.py`, `backend/scripts/build_similarity_index.py`, the `/api/similarity/*` routes in `backend/app.py`, and `frontend/src/tabs/SimilarityTab.tsx`. No content here was reconstructed from memory of past conversation — every claim, including the performance figures in §2.3, traces to a specific line of code or code comment read during this pass.
