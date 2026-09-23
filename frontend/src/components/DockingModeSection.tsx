import { useState } from "react";
import type { AdvancedDockingState } from "../lib/useAdvancedDocking";
import { SegmentedToggle } from "./SegmentedToggle";
import { BindingSiteModal } from "./BindingSiteModal";
import { TargetIcon } from "./Icons";

export function DockingModeSection({ adv, targetId }: { adv: AdvancedDockingState; targetId: string }) {
  const [modalOpen, setModalOpen] = useState(false);
  const { site, siteError, dockingMode, setDockingMode, selected, effectiveBox } = adv;

  let summary: React.ReactNode = <span className="text-inkmut">No binding-site evidence for this target.</span>;
  let showBtn = false;
  // Mode-first branching: a target can have blind_box_size available (any
  // receptor built on disk) while having NO site-specific center (the 5
  // targets neutralized to blind-only — see useAdvancedDocking's
  // loadBindingSite) — checking `site.center` here rather than mere
  // truthiness of `site` is what lets Blind mode show real numbers for
  // those targets instead of falling through to the "no site" branch
  // below just because a site-specific default doesn't exist.
  if (dockingMode === "blind") {
    if (site?.blind_box_size) {
      summary = (
        <>
          Blind docking: whole-protein search box{" "}
          <b>{site.blind_box_size.map((v) => v.toFixed(0)).join(" × ")}</b> Å — no pocket assumed.
        </>
      );
      showBtn = true;
    } else {
      summary = <span className="text-inkmut">Blind box unavailable — no prepared receptor on disk for this target.</span>;
    }
  } else if (site?.center) {
    // No separate "Automatic vs Manual" choice — the co-crystallized
    // ligand's pocket is selected by default the moment the site loads;
    // "View / edit binding site" below always opens the SAME view, listing
    // every residue in the receptor with the automatic pocket pre-checked,
    // so picking a different residue set or dragging the box are just
    // edits to that same default rather than a different starting mode.
    const [, activeSize] = effectiveBox();
    const total = site.allResidues.length || site.residues.length;
    summary = (
      <>
        Binding site: <b>{selected.size}</b> residue(s) selected (of {total} in the receptor) · box{" "}
        {(activeSize || site.box_size)?.map((v) => v.toFixed(1)).join(" × ")} Å
      </>
    );
    showBtn = true;
  } else if (targetId.startsWith("GENE_")) {
    summary = (
      <span className="text-inkmut">
        {siteError && /no validated small-molecule binding site/i.test(siteError)
          ? "No validated small-molecule binding site — no real co-crystallized ligand exists for this target."
          : "No automatic default yet for this target — pick a structure below (Advanced Settings)."}
      </span>
    );
  }

  return (
    <div>
      <label className="field-label">Docking mode</label>
      <SegmentedToggle
        value={dockingMode}
        onChange={(v) => setDockingMode(v as any)}
        options={[
          { value: "site_specific", label: "Site-specific" },
          { value: "blind", label: "Blind (whole protein)" },
        ]}
      />
      <div className="field-hint">{summary}</div>
      {showBtn && (
        <button type="button" className="btn-link mt-1.5" onClick={() => setModalOpen(true)}>
          <TargetIcon className="h-3.5 w-3.5" />
          View / edit binding site
        </button>
      )}
      {modalOpen && <BindingSiteModal adv={adv} targetId={targetId} onClose={() => setModalOpen(false)} />}
    </div>
  );
}
