import { useCallback, useEffect, useRef, useState } from "react";
import * as api from "./api";
import type {
  AdvancedDockingBody,
  BindingSiteResponse,
  PdbLigand,
  PocketResidue,
  ReceptorProfile,
  StructureCandidate,
} from "./types";

export type DockingMode = "site_specific" | "blind";

export interface SiteState {
  center?: [number, number, number];
  box_size?: [number, number, number];
  residues: PocketResidue[];
  /** The full receptor residue list (every amino acid), for manual
      binding-site selection — `residues` above stays the automatically-
      detected pocket neighborhood only. Falls back to `residues` itself
      when unavailable (an older cached profile, or a structure whose
      cleaned receptor is no longer on disk). */
  allResidues: PocketResidue[];
  blind_center?: [number, number, number];
  blind_box_size?: [number, number, number];
  receptorUrl: string | null;
  ligandUrl: string | null;
}

export const isGeneOnly = (targetId: string) => targetId.startsWith("GENE_");
export const residueKey = (r: PocketResidue) => `${r.chain}:${r.resnum}`;

/** Mirrors the original static/index.html's per-group `_advState[g]` blob —
    one instance per "Screen" / "Docking" tab, driven by that tab's current
    target selection. Owns binding-site evidence, the manual structure
    picker, residue-driven / drag-driven box overrides, and the
    exhaustiveness/poses/GNINA knobs — everything that feeds
    getAdvanced() -> AdvancedDocking on submit. */
export function useAdvancedDocking(targetId: string) {
  const [dockingMode, setDockingMode] = useState<DockingMode>("site_specific");
  const [site, setSite] = useState<SiteState | null>(null);
  const [siteError, setSiteError] = useState<string | null>(null);
  const [customProfile, setCustomProfile] = useState<ReceptorProfile | null>(null);
  const [selected, setSelected] = useState<Set<string>>(new Set());
  const [boxOverride, setBoxOverride] = useState<{ center: [number, number, number]; size: [number, number, number] } | null>(null);

  const [exhaustiveness, setExhaustiveness] = useState<string>("");
  const [nPoses, setNPoses] = useState<string>("");
  const [useGnina, setUseGnina] = useState(true);

  const [candidates, setCandidates] = useState<StructureCandidate[] | null>(null);
  const [candidatesLoading, setCandidatesLoading] = useState(false);
  const [pickedPdb, setPickedPdb] = useState<string | null>(null);
  const [structureStatus, setStructureStatus] = useState<{ kind: "muted" | "ok" | "warn" | "err"; text: string } | null>(null);
  const [preparingStructure, setPreparingStructure] = useState(false);

  // Co-crystallized ligands of whichever PDB is currently selected (via
  // selectPdb below) — NO ligand is pre-picked here; which one is
  // biologically relevant depends on the structure, not a generic
  // heuristic, so the user always chooses explicitly (see pickLigand).
  const [ligandOptions, setLigandOptions] = useState<PdbLigand[] | null>(null);
  const [ligandOptionsLoading, setLigandOptionsLoading] = useState(false);
  const [ligandOptionsError, setLigandOptionsError] = useState<string | null>(null);

  // What last produced the active box override, if any — lets a manual
  // drag WIN over (and stay put despite) a later residue-checkbox toggle;
  // without this, toggling any one residue after dragging would silently
  // recompute the box from the checked residues and discard the drag (the
  // "drag editing loses state" bug). A ref mirrors the state for the
  // synchronous read inside applyResidueSelection's debounced callback
  // (state alone can be one render stale there).
  const [boxSource, setBoxSourceState] = useState<"residues" | "drag" | null>(null);
  const boxSourceRef = useRef<"residues" | "drag" | null>(null);
  const setBoxSource = useCallback((v: "residues" | "drag" | null) => {
    boxSourceRef.current = v;
    setBoxSourceState(v);
  }, []);

  const residueTimer = useRef<ReturnType<typeof setTimeout> | null>(null);

  const applyResidueSelection = useCallback(
    async (sel: Set<string>, siteVal: SiteState | null, profile: ReceptorProfile | null) => {
      if (!siteVal) return;
      // A manual drag is a more specific, deliberate action than an
      // incidental residue toggle — once the user has dragged the box,
      // leave it alone until they explicitly reset (Reset to Automatic).
      if (boxSourceRef.current === "drag") return;
      const pool = siteVal.allResidues?.length ? siteVal.allResidues : siteVal.residues || [];
      const residues = pool.filter((r) => sel.has(residueKey(r)));
      if (!residues.length) {
        setBoxOverride(null);
        setBoxSource(null);
        return;
      }
      try {
        const body: any = {
          target_id: targetId,
          residues: residues.map((r) => ({ chain: r.chain, resnum: r.resnum })),
          padding: 8.0,
        };
        if (profile?.receptor_pdb) body.receptor_pdb = profile.receptor_pdb;
        const r = await api.boxFromResidues(body);
        setBoxOverride({ center: r.center, size: r.box_size });
        setBoxSource("residues");
      } catch {
        /* leave prior override in place — a transient failure here shouldn't
           disturb whatever box was already in effect */
      }
    },
    [targetId, setBoxSource]
  );

  /** Dragging a handle in the 3D view — always wins over whatever produced
      the box before (see boxSource above). Exposed to BindingSiteModal
      instead of the raw setBoxOverride setter so every drag is correctly
      tagged and can't be silently clobbered by a later residue toggle. */
  const setBoxOverrideFromDrag = useCallback(
    (center: [number, number, number], size: [number, number, number]) => {
      setBoxSource("drag");
      setBoxOverride({ center, size });
    },
    [setBoxSource]
  );

  const toggleResidue = useCallback(
    (key: string, checked: boolean) => {
      setSelected((prev) => {
        const next = new Set(prev);
        if (checked) next.add(key);
        else next.delete(key);
        if (residueTimer.current) clearTimeout(residueTimer.current);
        residueTimer.current = setTimeout(() => applyResidueSelection(next, site, customProfile), 400);
        return next;
      });
    },
    [applyResidueSelection, site, customProfile]
  );

  const loadStructureCandidates = useCallback(async (tid: string) => {
    setCandidatesLoading(true);
    setCandidates(null);
    setStructureStatus(null);
    setPickedPdb(null);
    try {
      const d = await api.structureCandidates(tid);
      setCandidates(d.candidates || []);
    } catch {
      setCandidates([]);
    } finally {
      setCandidatesLoading(false);
    }
  }, []);

  const loadBindingSite = useCallback(async (tid: string) => {
    setSelected(new Set());
    setBoxOverride(null);
    setBoxSource(null);
    setSiteError(null);
    if (!tid) {
      setSite(null);
      return;
    }
    try {
      const d: BindingSiteResponse = await api.bindingSite(tid);
      // No center means no validated site-specific default exists (e.g. a
      // target neutralized to blind-only because its only recorded
      // "ligand" turned out to be a buffer/phasing ion, see
      // scripts/audit_ligand_selection.py) — treat the same as a failed
      // call rather than showing a hollow "0 pocket residues" site.
      if (!d.center) {
        // No SITE-SPECIFIC default (e.g. one of the 5 targets neutralized
        // to blind-only because their only recorded "ligand" turned out to
        // be a buffer/phasing ion — see scripts/audit_ligand_selection.py),
        // but the backend computes blind_center/blind_box_size independent
        // of that (RP.box_from_receptor only needs the receptor file on
        // disk) — keep a SiteState around (center/box_size left undefined)
        // instead of discarding it via setSite(null), so Blind mode still
        // shows real numbers instead of falling through to the same "no
        // site" message shown for Site-specific mode.
        setSite({
          residues: [],
          allResidues: d.all_residues?.length ? d.all_residues : [],
          blind_center: d.blind_center,
          blind_box_size: d.blind_box_size,
          receptorUrl: api.apiUrl(`/api/docking/receptor/${tid}`),
          ligandUrl: null,
        });
        // Unlike the catch block below (a failed/404 lookup — no registry
        // entry at all, the normal state for most GENE_ targets), reaching
        // here means the profile itself loaded fine but has no site to
        // report — the backend always gives a specific reason for that
        // (e.g. "no validated small-molecule binding site... use Blind
        // mode" for the 5 neutralized targets), so show it regardless of
        // isGeneOnly; only the truly-generic fallback stays suppressed for
        // gene-only targets (most of which simply never had ANY structure).
        setSiteError(d.error || (isGeneOnly(tid) ? null : "No binding-site evidence for this target."));
        return;
      }
      const s: SiteState = {
        center: d.center,
        box_size: d.box_size,
        residues: d.pocket_residues || [],
        allResidues: d.all_residues?.length ? d.all_residues : d.pocket_residues || [],
        blind_center: d.blind_center,
        blind_box_size: d.blind_box_size,
        receptorUrl: api.apiUrl(`/api/docking/receptor/${tid}`),
        ligandUrl: d.has_reference_ligand_mol ? api.apiUrl(`/api/targets/${tid}/reference_ligand.sdf`) : null,
      };
      setSite(s);
      setSelected(new Set(s.residues.map(residueKey)));
    } catch {
      setSite(null);
      if (!isGeneOnly(tid)) setSiteError("No binding-site evidence for this target.");
    }
  }, []);

  // reset + reload whenever the owning tab's target selection changes
  useEffect(() => {
    setCustomProfile(null);
    setExhaustiveness("");
    setNPoses("");
    setUseGnina(true);
    setLigandOptions(null);
    setLigandOptionsError(null);
    if (!targetId) {
      setSite(null);
      setCandidates(null);
      return;
    }
    loadBindingSite(targetId);
    loadStructureCandidates(targetId);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [targetId]);

  /** Picking a structure (Advanced Settings' candidate list, or the
      "manual structure" override) — the resulting site (centered on THAT
      structure's own co-crystallized ligand, with all its residues listed
      and the pocket ones pre-checked) becomes the new default the same way
      an automatically-loaded target's site does; residue toggling / box
      dragging both remain available afterward exactly as for any target. */
  const applyBindingSiteFromProfile = useCallback((profile: ReceptorProfile) => {
    const s: SiteState = {
      center: profile.center,
      box_size: profile.box_size,
      residues: profile.binding_site_residues || [],
      allResidues: profile.all_residues?.length ? profile.all_residues : profile.binding_site_residues || [],
      blind_center: profile.blind_center,
      blind_box_size: profile.blind_box_size,
      receptorUrl: api.apiUrl(`/api/docking/receptor_file?path=${encodeURIComponent(profile.receptor_pdb)}`),
      ligandUrl: null,
    };
    setSite(s);
    setSelected(new Set(s.residues.map(residueKey)));
    setBoxOverride(null);
    setBoxSource(null);
  }, [setBoxSource]);

  /** Step 1 of the manual override: pick WHICH PDB to consider, and list
      every real co-crystallized ligand in it — nothing gets built yet.
      Picking a specific ligand (pickLigand below) is a separate, explicit
      step; there is no automatic "main ligand" applied here. */
  const selectPdb = useCallback(async (pdbId: string) => {
    setPickedPdb(pdbId);
    setCustomProfile(null);
    setBoxOverride(null);
    setSelected(new Set());
    setStructureStatus(null);
    setLigandOptions(null);
    setLigandOptionsError(null);
    setLigandOptionsLoading(true);
    try {
      const sub = await api.listPdbLigands(pdbId);
      while (true) {
        await api.sleep(1000);
        const j = await api.pollRetry(() => api.pdbLigandsJob(sub.job_id));
        if (j.status === "done") {
          setLigandOptions(j.ligands || []);
          break;
        }
        if (j.status === "error") {
          setLigandOptionsError(j.error || `Could not list ligands for ${pdbId}.`);
          break;
        }
      }
    } catch (e: any) {
      setLigandOptionsError(e.message || "Error");
    } finally {
      setLigandOptionsLoading(false);
    }
  }, []);

  /** Step 2: build a receptor centered on ONE explicitly-chosen ligand
      (resname=null/chain=null for a blind-only pick — a candidate with no
      real ligand at all). */
  const pickLigand = useCallback(
    async (pdbId: string, resname: string | null, chain: string | null) => {
      setPickedPdb(pdbId);
      setBoxOverride(null);
      setSelected(new Set());
      setStructureStatus({ kind: "muted", text: `Preparing receptor for ${pdbId}… (strip/repair/PDBQT, ~1–2 min)` });
      setPreparingStructure(true);
      try {
        const sub = await api.submitCustomReceptor({
          target_id: targetId,
          pdb_id: pdbId,
          chain: chain || undefined,
          ligand_resname: resname || undefined,
        });
        while (true) {
          await api.sleep(1500);
          const j = await api.pollRetry(() => api.customReceptorJob(sub.job_id));
          if (j.status === "running" && j.step) {
            setStructureStatus({ kind: "muted", text: `${j.step}…` });
          }
          if (j.status === "done" && j.profile) {
            setCustomProfile(j.profile);
            const p = j.profile;
            setStructureStatus({
              kind: "ok",
              text: resname ? `Using ${pdbId} — ${resname} (manual).` : `Using ${pdbId} (manual) — no ligand selected, Blind mode only.`,
            });
            applyBindingSiteFromProfile(p);
            break;
          }
          if (j.status === "error") {
            setCustomProfile(null);
            setStructureStatus({ kind: "err", text: `Could not prepare ${pdbId}: ${j.error || "unknown error"}` });
            break;
          }
        }
      } catch (e: any) {
        setCustomProfile(null);
        setStructureStatus({ kind: "err", text: e.message || "Error" });
      } finally {
        setPreparingStructure(false);
      }
    },
    [targetId, applyBindingSiteFromProfile]
  );

  const resetToAutomatic = useCallback(() => {
    setExhaustiveness("");
    setNPoses("");
    setUseGnina(true);
    setCustomProfile(null);
    setLigandOptions(null);
    setLigandOptionsError(null);
    setPickedPdb(null);
    setBoxOverride(null);
    setBoxSource(null);
    setStructureStatus(null);
    if (targetId) loadBindingSite(targetId);
  }, [targetId, loadBindingSite, setBoxSource]);

  /** Whole-protein box in blind mode (that toggle IS the explicit choice —
      there's no pocket to define either way); else the box override (drag
      or residue selection) if present; else the target's automatic
      ligand-centered default the instant a site loads — no separate
      confirmation step, "automatic" is simply what's active until the
      user changes it via residue toggles or dragging. */
  const effectiveBox = useCallback((): [[number, number, number] | undefined, [number, number, number] | undefined] => {
    if (!site) return [undefined, undefined];
    if (dockingMode === "blind") return [site.blind_center || site.center, site.blind_box_size || site.box_size];
    if (boxOverride) return [boxOverride.center, boxOverride.size];
    return [site.center, site.box_size];
  }, [site, dockingMode, boxOverride]);

  const getAdvanced = useCallback((): AdvancedDockingBody | null => {
    const adv: AdvancedDockingBody = {};
    const exh = parseFloat(exhaustiveness);
    if (exhaustiveness !== "" && !isNaN(exh)) adv.exhaustiveness = Math.round(exh);
    const poses = parseFloat(nPoses);
    if (nPoses !== "" && !isNaN(poses)) adv.n_poses = Math.round(poses);
    if (!useGnina) adv.use_gnina = false;
    if (boxOverride) {
      adv.box_center = boxOverride.center;
      adv.box_size = boxOverride.size;
    }
    if (customProfile) adv.custom_profile = customProfile;
    if (dockingMode === "blind") adv.docking_mode = "blind";
    return Object.keys(adv).length ? adv : null;
  }, [exhaustiveness, nPoses, useGnina, boxOverride, customProfile, dockingMode]);

  return {
    targetId,
    dockingMode,
    setDockingMode,
    site,
    siteError,
    customProfile,
    selected,
    toggleResidue,
    boxOverride,
    boxSource,
    setBoxOverrideFromDrag,
    exhaustiveness,
    setExhaustiveness,
    nPoses,
    setNPoses,
    useGnina,
    setUseGnina,
    candidates,
    candidatesLoading,
    pickedPdb,
    structureStatus,
    preparingStructure,
    selectPdb,
    pickLigand,
    ligandOptions,
    ligandOptionsLoading,
    ligandOptionsError,
    resetToAutomatic,
    effectiveBox,
    getAdvanced,
  };
}

export type AdvancedDockingState = ReturnType<typeof useAdvancedDocking>;
