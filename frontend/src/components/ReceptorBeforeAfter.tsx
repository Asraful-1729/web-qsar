import { ReceptorPreview } from "./ReceptorPreview";
import { apiUrl } from "../lib/api";

/** Side-by-side original-vs-prepared receptor comparison for a manually
    picked Advanced Settings structure — the "before/after 3D preview" part
    of receptor-prep transparency (the step-by-step pipeline labels are
    surfaced live during prep itself, see useAdvancedDocking's pickLigand).
    Reuses ReceptorPreview (style switcher included) for both sides rather
    than a bespoke viewer. */
export function ReceptorBeforeAfter({
  rawPdbPath,
  cleanPdbPath,
  heightClassName = "h-[240px]",
}: {
  rawPdbPath: string;
  cleanPdbPath: string;
  /** Tailwind height class for each side's viewer canvas — a small fixed
      height inline, or e.g. "h-full" to fill a modal. */
  heightClassName?: string;
}) {
  const rawUrl = apiUrl(`/api/docking/receptor_file?path=${encodeURIComponent(rawPdbPath)}`);
  const cleanUrl = apiUrl(`/api/docking/receptor_file?path=${encodeURIComponent(cleanPdbPath)}`);
  return (
    <div className="grid h-full grid-cols-1 gap-2.5 sm:grid-cols-2">
      <div className="flex h-full flex-col">
        <div className="field-hint mb-1 shrink-0">Before — original PDB (waters, heteroatoms, missing atoms as deposited)</div>
        <div className="min-h-0 flex-1">
          <ReceptorPreview receptorUrl={rawUrl} ligandUrl={null} heightClassName={heightClassName} showHetero />
        </div>
      </div>
      <div className="flex h-full flex-col">
        <div className="field-hint mb-1 shrink-0">After — stripped, repaired, protonated (this is what gets docked against)</div>
        <div className="min-h-0 flex-1">
          <ReceptorPreview receptorUrl={cleanUrl} ligandUrl={null} heightClassName={heightClassName} />
        </div>
      </div>
    </div>
  );
}
