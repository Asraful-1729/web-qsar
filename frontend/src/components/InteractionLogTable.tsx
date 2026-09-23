import { interactionLogRows, interactionLogCsv } from "../lib/interactionLog";
import type { DockResultRow } from "../lib/types";

function downloadCsv(csv: string, name: string) {
  const url = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  URL.revokeObjectURL(url);
}

type Results = (DockResultRow | null | undefined)[];

/** The collapsible toggle for InteractionLogTable below, split out so it can
    sit inline in the results toolbar (left of Export) instead of floating
    as its own row — open/close state is owned by the caller so both pieces
    stay in sync despite living in different parts of the layout. */
export function InteractionTableToggle({ results, open, onToggle }: { results: Results; open: boolean; onToggle: () => void }) {
  const rows = interactionLogRows(results);
  if (!rows.length) return null;
  return (
    <button type="button" className="btn-link" onClick={onToggle}>
      {open ? "▾" : "▸"} Interaction table ({rows.length} interaction{rows.length === 1 ? "" : "s"})
    </button>
  );
}

/** Cross-compound interaction log: one row per interaction EVENT (never
    aggregated by residue or compound), so it's immediately clear which
    compound is interacting with which residue, through what bond type,
    and at what distance. Rendered only when `open` (see
    InteractionTableToggle) — a secondary view next to the main
    per-compound results table, not the primary one. */
export function InteractionLogTable({ results, fileBaseName, open }: { results: Results; fileBaseName: string; open: boolean }) {
  const rows = interactionLogRows(results);
  if (!open || !rows.length) return null;
  return (
    <div className="border-t border-line px-5 py-3.5">
      <div className="mb-2 flex justify-end">
        <button type="button" className="btn-link" onClick={() => downloadCsv(interactionLogCsv(rows), `${fileBaseName}_interactions.csv`)}>
          Download CSV
        </button>
      </div>
      <div className="max-h-[320px] overflow-y-auto rounded-lg border border-line">
        <table className="w-full border-collapse text-[12.5px]">
          <thead>
            <tr>
              {["Compound", "Interacting Residue", "Bond Type", "Length"].map((h) => (
                <th key={h} className="sticky top-0 z-10 border-b border-line bg-surface2 px-2.5 py-2 text-left text-[10.5px] font-semibold uppercase tracking-wide text-inkmut">
                  {h}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((r, i) => (
              <tr key={i} className="hover:bg-canvas">
                <td className="smi-mono max-w-[240px] border-b border-surface2 px-2.5 py-1.5 text-ink" title={r.compound}>{r.compound}</td>
                <td className="border-b border-surface2 px-2.5 py-1.5 font-semibold text-ink">{r.residue}</td>
                <td className="border-b border-surface2 px-2.5 py-1.5 text-inkmut">{r.bondType}</td>
                <td className="border-b border-surface2 px-2.5 py-1.5 text-inkmut">{r.length != null ? `${r.length} Å` : "—"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
