import { useEffect, useState } from "react";
import * as api from "../lib/api";
import type { RecommendationResponse } from "../lib/types";

export function WhyThisButton({ targetId }: { targetId: string }) {
  const [rec, setRec] = useState<RecommendationResponse | "loading" | "none" | null>(null);

  // Drop any previous target's evidence the moment the target changes,
  // rather than relying on a remount-via-key (React 18 has been observed
  // to leave the OLD keyed instance mounted alongside the new one here —
  // a real duplicate-DOM bug, not just stale state — so this component is
  // never remounted on target change; it just resets itself instead).
  useEffect(() => {
    setRec(null);
  }, [targetId]);

  const show = async () => {
    setRec("loading");
    try {
      setRec(await api.targetRecommendation(targetId));
    } catch {
      setRec("none");
    }
  };

  return (
    <div>
      <button type="button" className="btn-ghost mt-2" onClick={show}>
        Why this?
      </button>
      {rec === "loading" && <div className="field-hint">Loading…</div>}
      {rec === "none" && <div className="field-hint">No structural-evidence data for this target.</div>}
      {rec && typeof rec === "object" && <RecommendationCard rec={rec} />}
    </div>
  );
}

export function RecommendationCard({ rec }: { rec: RecommendationResponse }) {
  const s = rec.structure;
  const p = rec.panel_evidence;
  return (
    <div className="mt-2 rounded-xl border border-line bg-surface2/40 px-3 py-2.5">
      <div className="mb-1 text-[13px]">{rec.headline}</div>
      {s.pdb_id && (
        <div className="text-[12.5px] text-ink/80">
          PDB {s.pdb_id}/{s.ligand_resname || ""}
          {s.rank_in_panel_evidence != null ? <> &nbsp;·&nbsp; rank #{s.rank_in_panel_evidence} in panel evidence</> : null}
        </div>
      )}
      <h5 className="mb-1 mt-2.5 text-[11px] font-bold uppercase tracking-wide text-brand-700">Structural evidence (panel_results_v2.csv)</h5>
      <div className="text-[12.5px] text-ink/80">
        Top-ranked: <b>{p.top_ranked_pdb_id || "—"}{p.top_ranked_chain ? ":" + p.top_ranked_chain : ""}</b> ({p.top_ranked_ligand || "?"})
        <br />
        Resolution {p.resolution ?? "—"} Å ({p.resolution_tier || "?"}) · RSCC {p.ligand_RSCC ?? "—"} · RSR {p.ligand_RSR ?? "—"} · R-free {p.r_free ?? "—"}
        <br />
        {p.n_qualifying_structures ?? "?"} qualifying structures · {p.chembl_activity_records ?? "?"} ChEMBL activity records
      </div>
      {p.note && <div className="mt-1.5 text-[11px] text-inkmut">{p.note}</div>}
    </div>
  );
}
