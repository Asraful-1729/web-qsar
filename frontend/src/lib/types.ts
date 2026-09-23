// Loose but useful TypeScript shapes mirroring app.py's JSON responses.
// Deeply-variable nested payloads (docking results, ADMET learned groups,
// compare views) are typed permissively — the UI reads through them
// defensively exactly like the original vanilla-JS client did.

export interface TargetMeta {
  target_id: string;
  name?: string;
  best_model?: string | null;
  n_compounds?: number | null;
  test_r2?: number | null;
  test_rmse?: number | null;
  ad_coverage_pct?: number | null;
  tropsha_pass?: boolean | null;
}

export interface DiseaseSummary {
  disease_id: string;
  name: string;
  is_therapeutic_area?: boolean;
}

export interface DiseaseTarget {
  target_id?: string;
  target_symbol: string;
  has_qsar_model: boolean;
  /** Whether the model bucket is actually downloaded on this machine right
      now — has_qsar_model is download-independent (true for every QSAR-
      panel gene); this tells the picker whether selecting it will need to
      fetch it first (see useDownloadGate.ts). */
  model_installed?: boolean;
  disease_score: number | string;
}

export type Confidence = "high" | "med" | "low" | "out" | "na";

export interface QsarRow {
  input_smiles?: string;
  smiles?: string;
  parsed_ok?: boolean;
  predicted_pIC50?: number | null;
  in_domain: boolean;
  ad_z?: number | null;
  confidence?: Confidence;
  confidence_label?: string;
  confidence_basis?: string;
  rank?: number;
}

export interface PredictResponse {
  target: { id: string; name: string };
  model?: string | null;
  model_metrics: {
    test_r2?: number | null;
    test_rmse?: number | null;
    pearson_r?: number | null;
    ad_coverage_pct?: number | null;
    tropsha_pass?: boolean | null;
    y_random_delta_r2?: number | null;
  };
  counts: { in_domain: number; out_of_domain: number; skipped: number; submitted: number };
  in_domain: QsarRow[];
  out_of_domain: QsarRow[];
  skipped: string[];
  disclaimer: string;
}

export interface AdmetProfile {
  input_smiles: string;
  standardised_smiles?: string;
  parsed_ok: boolean;
  physicochemical?: { mw: number; logp: number; qed: number };
  drug_likeness_flags?: { lipinski_pass: boolean; lipinski_violations: number };
  n_alerts?: number;
  learned?: {
    available: boolean;
    note?: string;
    flags?: any[];
    groups?: Record<
      string,
      { name: string; label: string; task: "class" | "reg"; unit: string; display: string; tone: string; percentile?: number | null }[]
    >;
    /** Every column ADMET-AI actually returned for this compound (all
        physicochemical descriptors, task predictions, and
        *_drugbank_approved_percentile fields) — groups/flags above are
        only the curated subset shown in the UI. Used for the "Download
        all data (CSV)" export. */
    raw?: Record<string, number | string | boolean | null>;
  };
}

export interface AdmetResponse {
  mode?: "result" | "job";
  job_id?: string;
  total?: number;
  status?: string;
  done?: number;
  profiles: AdmetProfile[];
  learned: { available: boolean; note?: string };
  disclaimer: string;
}

export interface CompareResponse {
  disclaimer: string;
  targets: { target_id: string }[];
  matrix: { smiles?: string; input_smiles?: string; cells: Record<string, { pred: number | null; in_domain: boolean }>; coverage: { in_domain_targets: number; total_targets: number } }[];
  consensus_ranking: { consensus_rank: number; smiles: string; n_active: number; mean_pred: number; coverage: string }[];
  selective_candidates: { smiles: string; target: string; pred: number; gap: number }[];
  selective_gap: number;
  multi_target_candidates: { smiles: string; active_targets: string[]; n_active: number; mean_pred: number }[];
  active_cut: number;
  best_per_target: Record<string, { smiles: string; pred: number }[]>;
  admet: Record<string, AdmetProfile>;
}

export interface AdvancedDockingBody {
  exhaustiveness?: number | null;
  n_poses?: number | null;
  use_gnina?: boolean | null;
  docking_mode?: "site_specific" | "blind" | null;
  box_center?: [number, number, number] | null;
  box_size?: [number, number, number] | null;
  custom_profile?: ReceptorProfile | null;
}

export interface ReceptorProfile {
  target_id?: string;
  receptor_pdb: string;
  center: [number, number, number];
  box_size: [number, number, number];
  pdb_source?: string;
  /** Original, unmodified structure as fetched from RCSB — before
      strip/repair. Lets the UI show a before/after comparison against
      receptor_pdb (the cleaned one). Only set for a manual Advanced
      Settings structure pick, not the registry's build-time default. */
  raw_pdb_path?: string;
  binding_site_residues?: PocketResidue[];
  /** Full receptor residue list — see BindingSiteResponse's field of the
      same name. */
  all_residues?: PocketResidue[];
  blind_center?: [number, number, number];
  blind_box_size?: [number, number, number];
  /** B13 — real before/after facts for each receptor-prep stage (atom
      counts, pocket residues found, box dimensions) — only set for a
      manual Advanced Settings structure pick, same as raw_pdb_path,
      since the registry's pre-built defaults don't carry this record. */
  prep_report?: { label: string; detail: string }[];
  /** Real wall-clock seconds per phase, so a slow "prepare structure" run
      is diagnosable instead of one opaque timed-but-unbroken-down wait. */
  timing?: {
    fetch_pdb_seconds?: number;
    build_seconds?: number;
  };
  reference_ligand_resname?: string;
  /** Set only when the true PDB Chemical Component Dictionary code (e.g.
      "A1AWR") had to be truncated to fit into reference_ligand_resname's
      legacy-PDB-compatible 3-character form (e.g. "A1A") — see
      pdb_fetch.py's resolve_ligand_resname. Purely informational/display;
      lookups must keep using reference_ligand_resname. */
  reference_ligand_full_id?: string | null;
  /** Heavy-atom-only MW estimate (Da) for the selected reference ligand —
      a coarse proxy, not a cheminformatics-grade value; see
      docking/receptor_prep.py's _approx_heavy_atom_mw. */
  reference_ligand_mw?: number;
  /** Non-blocking flags from receptor_prep.py's ligand sanity checks
      (too small/large, too few pocket contacts) — surfaced so a
      suspicious pick is visible, never silently accepted or rejected. */
  ligand_sanity_warnings?: string[];
  /** "crystallographic_annotation" (a specific ligand was named — curated
      data or a user's own pick) vs. "automatic" (nobody knew, fell back to
      "largest non-additive HETATM group"). */
  ligand_selection_method?: "crystallographic_annotation" | "automatic";
  /** "ok" vs. "review_required" (ligand_sanity_warnings non-empty) —
      "buildable" is not the same claim as "biologically validated". */
  build_status?: "ok" | "review_required";
  [k: string]: any;
}

/** One real co-crystallized ligand found in a job's raw PDB structure —
    see backend/docking/receptor_prep.py's list_ligands(). */
export interface AlternateLigand {
  resname: string;
  chain: string;
  resnum: number;
  n_atoms: number;
  center: [number, number, number];
}
export interface AlternateLigandsResponse {
  available: boolean;
  note?: string;
  current: AlternateLigand | null;
  ligands: AlternateLigand[];
}

export interface EnrichmentRunSettings {
  pdb_source?: string | null;
  decoy_method?: string | null;
  engine?: string | null;
  exhaustiveness?: number | null;
  center?: [number, number, number] | null;
  box_size?: [number, number, number] | null;
  docking_mode?: "site_specific" | "blind" | null;
}
export interface PocketResidue {
  chain: string;
  resnum: number;
  resname: string;
}

export interface BindingSiteResponse {
  target_id: string;
  center?: [number, number, number];
  box_size?: [number, number, number];
  reference_ligand_resname?: string;
  pocket_residues: PocketResidue[];
  n_pocket_residues: number;
  /** The FULL receptor residue list (every amino acid, not just the
      auto-detected pocket neighborhood) — what a manual binding-site
      picker selects from. */
  all_residues?: PocketResidue[];
  has_reference_ligand_mol: boolean;
  blind_center?: [number, number, number];
  blind_box_size?: [number, number, number];
  error?: string | null;
}

/** Redocking-pose validation: re-dock the target's own known reference
    ligand and compare the top pose to its experimentally-observed (crystal)
    pose via atom-map-safe RMSD — see docking/pipeline.py's
    redock_reference_for_profile. Runs AUTOMATICALLY as part of a docking/
    screen submission (same PDB/exhaustiveness/n_poses the user picked for
    their compounds) — there is no separate action to trigger it, and it
    never gates or blocks the compound-docking results it rides along with. */
export interface RedockingValidationResult {
  status: "ok" | "engine_unavailable" | "ligand_prep_failed" | "no_pose" | "no_crystal_reference" | "error";
  engine?: string;
  top_score?: number;
  n_poses?: number;
  redocked_pose_pdb?: string;
  /** The experimental (crystal) pose, verbatim SDF — same molecule as
      redocked_pose_pdb, for the "experimental vs redocked" 3D overlay.
      Only present alongside a real reference_rmsd (both come from the
      same crystal_sdf file). */
  crystal_pose_sdf?: string;
  reference_rmsd?: number;
  validated?: boolean;
  rmsd_error?: string;
  reference_ligand_resname?: string;
  error?: string;
}

/** One real co-crystallized ligand found in a PDB entry — see
    receptor_prep.list_ligands. Deliberately carries no "is this the main
    one" flag: which ligand matters depends on the structure's own biology,
    not a generic heuristic, so the picker shows every candidate and lets
    the user choose explicitly. */
export interface PdbLigand {
  resname: string;
  chain: string;
  resnum: number;
  n_atoms: number;
  center: [number, number, number];
}

export interface PdbLigandsJobStatus {
  status: "queued" | "running" | "done" | "error";
  step?: string;
  ligands?: PdbLigand[];
  error?: string;
}

export interface StructureCandidate {
  pdb_id: string;
  /** null only for a blind_only fallback candidate (see below) — every
      gene's own real ligand candidates always have a real resname. */
  resname: string | null;
  csv_rank: number | null;
  is_current_default: boolean;
  resolution?: number | null;
  ligand_RSCC?: number | null;
  ligand_RSR?: number | null;
  /** Set when this gene's own real ligand candidates were all filtered
      out (buffer/cryoprotectant/ion, or below the MW floor — see
      app.py's _structure_candidates_for_gene) and this is the best-
      quality structure on record instead, offered for Blind docking
      only — picking it builds a receptor with no site-specific default
      (docking/receptor_prep.py's build_receptor() no-ligand fallback). */
  blind_only?: boolean;
}

export interface StructureCandidatesResponse {
  target_id?: string;
  gene: string;
  default_pdb_id: string | null;
  candidates: StructureCandidate[];
  n_qualifying_structures?: number | null;
  note?: string;
}

export interface DockResultRow {
  smiles: string;
  status: string;
  reason?: string;
  error?: string;
  /** Finer failure classification than `status` (e.g. "invalid_molecule",
      "conformer_generation_failed", "pose_validity_failed") plus a plain-
      language suggested_action — absent for a successful ("ok") result,
      since there's nothing to diagnose. See backend/docking/
      failure_diagnostics.py for the exact category list. */
  category?: string;
  suggested_action?: string;
  confidence?: "high" | "medium" | "low" | "none";
  vina_score?: number | null;
  n_valid?: number;
  pose_self_consistency?: number;
  gnina?: { cnn_score?: number | null; cnn_affinity?: number | null; gnina_affinity?: number | null };
  interaction_png?: string | null;
  interaction_source?: string;
  residue_overlap_pct?: number | null;
  interactions?: {
    name?: string;
    residue?: string;
    resname?: string;
    resid?: number;
    chain?: string;
    category?: string;
    label?: string;
    type?: string;
    distance?: number;
  }[];
  pose_pdb?: string | null;
}

export interface DockJobDone {
  status: "done" | "cancelled";
  done: number;
  total: number;
  caveat?: string | null;
  results: DockResultRow[];
  receptor_pdb_path?: string | null;
  pdb_source?: string | null;
  redocking_validation?: RedockingValidationResult | null;
}
export interface DockJobPending {
  status: "queued" | "running";
  done: number;
  total: number;
  caveat?: string | null;
  pdb_source?: string | null;
}
export interface DockJobError {
  status: "error";
  error?: string;
}
export type DockJobStatus = DockJobDone | DockJobPending | DockJobError;

export interface ScreenShortlistRow {
  rank: number;
  input_smiles: string;
  smiles: string;
  qsar: QsarRow & { in_domain: boolean };
  docking: DockResultRow | null;
  fused_score?: number | null;
  caveats: string[];
}

export interface ScreenResult {
  target_id: string;
  counts: { submitted: number; parsed: number; skipped: number };
  methods_note: string;
  docking_used: boolean;
  docking_note?: string;
  shortlist: ScreenShortlistRow[];
  skipped: string[];
  receptor_pdb_path?: string | null;
  pdb_source?: string | null;
  redocking_validation?: RedockingValidationResult | null;
}

export interface ScreenJobStatus {
  status: "queued" | "running" | "done" | "error" | "cancelled";
  step?: number;
  step_label?: string;
  total_steps?: number;
  done?: number | null;
  total?: number | null;
  result?: ScreenResult;
  error?: string;
}

export interface RecommendationResponse {
  headline: string;
  structure: {
    pdb_id?: string;
    ligand_resname?: string;
    rank_in_panel_evidence?: number | null;
  };
  panel_evidence: {
    top_ranked_pdb_id?: string;
    top_ranked_chain?: string;
    top_ranked_ligand?: string;
    resolution?: number | null;
    resolution_tier?: string;
    ligand_RSCC?: number | null;
    ligand_RSR?: number | null;
    r_free?: number | null;
    n_qualifying_structures?: number | null;
    chembl_activity_records?: number | null;
    note?: string;
  };
}

export interface DockingStatus {
  ready: boolean;
  note?: string;
  import_error?: string;
  packages?: Record<string, boolean>;
  package_desc?: Record<string, string>;
  binaries?: Record<string, boolean>;
  binary_desc?: Record<string, string>;
  planned?: string[];
  docking_targets?: string[];
  target_details?: {
    target_id: string;
    name: string;
    site_source?: string;
  }[];
}

export interface BucketFile {
  path: string;
  name: string;
  bytes: number;
  annotation?: string;
  category: string;
}

// ---------- on-demand downloads (Downloads tab) ----------
export interface DownloadKindStatus {
  available: boolean;
  installed: boolean;
  size?: number;
}
export interface DownloadTargetRow {
  target_id: string;
  model: DownloadKindStatus;
  docking: DownloadKindStatus;
}
export interface DownloadsStatus {
  download_base_url: string | null;
  targets: DownloadTargetRow[];
}
export interface DownloadStartResponse {
  job_id: string | null;
  already_installed?: boolean;
}
export interface DownloadJobStatus {
  target_id: string;
  kind: "model" | "docking";
  state: "starting" | "downloading" | "extracting" | "done" | "error" | "cancelled";
  done: number;
  total: number;
  error?: string | null;
}

export interface SimilarityHit {
  id: string;
  smiles: string;
  name?: string | null;
  tanimoto: number;
  scaffold_match: boolean;
  molecular_weight?: number | null;
  chemical_super_class?: string | null;
  np_classifier_class?: string | null;
  mcs_smarts?: string;
  mcs_n_atoms?: number;
}
export interface SimilarityResult {
  query_scaffold?: string | null;
  n_indexed: number;
  n_results: number;
  results: SimilarityHit[];
}

export interface TargetFishingCompound {
  smiles: string;
  tanimoto: number;
  pchembl_value?: number | null;
}
export interface TargetFishingHit {
  target_chembl: string;
  /** Human-readable ChEMBL target name — always present; use this for
      display, target_id is not. */
  target_pref_name: string | null;
  /** Only set when this ChEMBL target happens to also be one of this
      app's own docking/QSAR targets — most target-fishing hits won't
      have one, since the reference pool is much broader than this app's
      own target list. null means "not one of our targets," not an error. */
  target_id: string | null;
  n_similar_actives: number;
  best_similarity: number;
  mean_similarity: number;
  best_pchembl: number | null;
  mean_pchembl: number | null;
  n_scaffolds: number;
  /** A documented, rule-based combination of similarity + scaffold
      diversity + potency (see backend/target_fishing.py's module
      docstring for the exact formula and weights) — NOT a probability,
      NOT statistically calibrated. Use for ranking/sorting only. */
  evidence_score: number;
  compounds: TargetFishingCompound[];
}
export interface TargetFishingResult {
  n_indexed_compound_target_pairs: number;
  n_targets_searched: number;
  n_targets_matched: number;
  results: TargetFishingHit[];
}
/** Target Prediction v2 — a separate, independently-validated method from
    TargetFishing (v1) above, not a drop-in replacement. See
    target_prediction_v2/METHODS_AND_VALIDATION.md for the full validation
    account (significantly beats v1 and a real SEA reimplementation on
    genuinely held-out data) and backend/target_prediction_v2.py's module
    docstring for the method itself. */
export interface TargetPredictionV2NativeNeighbour {
  smiles: string;
  tanimoto: number;
  pchembl_value: number | null;
  potency_weight: number;
}
export interface TargetPredictionV2OrthologueNeighbour {
  smiles: string;
  tanimoto: number;
  species_provenance: string | null;
}
export interface TargetPredictionV2Evidence {
  /** Present only in the "pooled" regime. */
  native_neighbours?: TargetPredictionV2NativeNeighbour[];
  orthologue_neighbours?: TargetPredictionV2OrthologueNeighbour[];
  /** Present only in the "best_similarity_fallback" regime. */
  best_compound_smiles?: string | null;
  best_compound_pchembl?: number | null;
}
export interface TargetPredictionV2Hit {
  rank: number;
  target_chembl: string;
  target_pref_name: string | null;
  /** L-score (pooled regime) or best_similarity (fallback regime) — NEVER
      a calibrated probability, see confidence_label for the exact caveat
      (Phase 4 measured calibration and it failed to beat this raw score). */
  score: number;
  confidence_label: string;
  evidence: TargetPredictionV2Evidence;
}
export interface TargetPredictionV2Result {
  /** Which retrieval regime fired for this query — density-adaptive gate,
      re-fit and G2-validated against the rebuilt v2 index. */
  regime: "pooled" | "best_similarity_fallback";
  density: number;
  density_threshold: number;
  n_neighbours_considered: number;
  n_targets_indexed: number;
  results: TargetPredictionV2Hit[];
}

export interface CuratedCompoundSuggestion {
  smiles: string;
  target_chembl: string;
  target_pref_name: string | null;
  target_id: string | null;
  pchembl_value?: number | null;
}

export interface LiteraturePaper {
  pmid: string | null;
  title: string;
  abstract?: string | null;
  journal?: string | null;
  year?: string | null;
  authors?: string;
  url?: string | null;
}
/** A5's evidence-chain report — deeply variable per-section shape
    (each section is independently "available: boolean" depending on
    what data existed for this compound), so typed permissively like
    the rest of this file's larger nested payloads. */
export interface ResearchReport {
  generated_at: string;
  target_id: string;
  pipeline_stages: string[];
  natural_source: { plant_source?: string | null };
  chemical_identity: Record<string, any>;
  reported_activity: { available: boolean; query?: string; n_results?: number; papers?: LiteraturePaper[]; note?: string };
  target_prediction: Record<string, any>;
  qsar_prediction: Record<string, any>;
  docking: DockResultRow | Record<string, any>;
  interaction_analysis: { interactions: any[]; residue_overlap_pct?: number | null };
  admet: { available: boolean; profile?: AdmetProfile; note?: string };
  off_target_analysis: Record<string, any>;
  evidence_summary: string;
  methods_draft: string;
  reproducibility: Record<string, any>;
  disclaimer: string;
}
export interface ResearchReportResponse {
  report: ResearchReport;
  markdown: string;
}
