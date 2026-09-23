import { useEffect, useState } from "react";
import type { AdvancedDockingState } from "../lib/useAdvancedDocking";
import { RefreshIcon } from "./Icons";

export function AdvancedSettingsPanel({
  adv,
  openByDefault = false,
}: {
  adv: AdvancedDockingState;
  openByDefault?: boolean;
}) {
  const [open, setOpen] = useState(openByDefault);

  // Reacts to openByDefault flipping true (a target just got picked)
  // instead of remounting via a `key` prop — see WhyThisButton's comment
  // for why: this component is never remounted on target change.
  useEffect(() => {
    if (openByDefault) setOpen(true);
  }, [openByDefault]);

  return (
    <div className="mt-3.5 rounded-xl border border-line">
      <button
        type="button"
        onClick={() => setOpen((o) => !o)}
        className="flex w-full items-center justify-between px-3 py-2.5 text-left text-[12.5px] font-semibold text-brand-700"
      >
        Advanced Settings
        <svg
          viewBox="0 0 24 24"
          width="14"
          height="14"
          fill="none"
          stroke="currentColor"
          strokeWidth="2.4"
          className={`transition-transform ${open ? "rotate-180" : ""}`}
        >
          <path d="M6 9l6 6 6-6" strokeLinecap="round" strokeLinejoin="round" />
        </svg>
      </button>
      {open && (
        <div className="border-t border-line px-3 pb-3 pt-2.5">
          <div className="mb-3">
            <label className="field-label">Exhaustiveness</label>
            <input
              type="number"
              min={1}
              max={64}
              placeholder="Automatic (8)"
              className="field-input"
              value={adv.exhaustiveness}
              onChange={(e) => adv.setExhaustiveness(e.target.value)}
            />
          </div>
          <div className="mb-3">
            <label className="field-label">Number of poses</label>
            <input
              type="number"
              min={1}
              max={20}
              placeholder="Automatic (9)"
              className="field-input"
              value={adv.nPoses}
              onChange={(e) => adv.setNPoses(e.target.value)}
            />
          </div>
          <div className="mb-3">
            <label className="flex cursor-pointer items-center gap-2 text-[12.5px] font-medium text-ink">
              <input type="checkbox" checked={adv.useGnina} onChange={(e) => adv.setUseGnina(e.target.checked)} />
              GNINA CNN rescoring (second opinion, if installed)
            </label>
          </div>
          <div className="mb-3">
            <label className="field-label">Binding box</label>
            <div className="field-hint">
              {adv.boxSource === "drag"
                ? "Manually adjusted — dragged in \"View / edit binding site\" above."
                : adv.boxSource === "residues"
                ? "Set from the residues selected in \"View / edit binding site\" above."
                : adv.site?.center
                ? "Automatic — centered on the co-crystallized reference ligand. Open \"View / edit binding site\" above to change residues or drag to adjust."
                : "Not yet defined — no validated binding site for this target."}
            </div>
          </div>
          <button type="button" className="btn-ghost" onClick={() => adv.resetToAutomatic()}>
            <RefreshIcon className="h-3.5 w-3.5" />
            Reset to Automatic
          </button>
        </div>
      )}
    </div>
  );
}
