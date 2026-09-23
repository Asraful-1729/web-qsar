import type { DockResultRow } from "./types";

export interface InteractionLogRow {
  /** The actual compound identifier (SMILES) — never a generic per-row
      label — so it's unambiguous which compound each interaction belongs
      to even when several compounds hit the same residue. */
  compound: string;
  compoundIdx: number;
  residue: string; // e.g. "PHE182" — resname+resid, or chain:resid if resname is missing
  bondType: string; // e.g. "Conventional Hydrogen Bond", "Pi-Alkyl", "Salt Bridge"
  length?: number; // Å
}

/** One row per interaction EVENT, never aggregated across compounds or
    interaction types — if PHE182 interacts with 3 different compounds,
    that's 3 separate rows, each naming its own compound. Replaces the
    former per-residue frequency summary (which hid which compound was
    responsible for which contact). */
export function interactionLogRows(results: (DockResultRow | null | undefined)[]): InteractionLogRow[] {
  const rows: InteractionLogRow[] = [];
  results.forEach((r, compoundIdx) => {
    if (!r?.interactions?.length) return;
    for (const h of r.interactions) {
      const residue = h.residue || (h.chain && h.resid != null ? `${h.chain}:${h.resid}` : "—");
      rows.push({
        compound: r.smiles,
        compoundIdx,
        residue,
        bondType: h.label || h.category || h.type || "—",
        length: h.distance,
      });
    }
  });
  return rows.sort((a, b) => a.compoundIdx - b.compoundIdx || (a.length ?? 0) - (b.length ?? 0));
}

function csvCell(v: unknown): string {
  if (v === null || v === undefined) return "";
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

export function interactionLogCsv(rows: InteractionLogRow[]): string {
  const header = ["compound", "interacting_residue", "bond_type", "length_angstrom"];
  const lines = [header.map(csvCell).join(",")];
  for (const r of rows) {
    lines.push([r.compound, r.residue, r.bondType, r.length].map(csvCell).join(","));
  }
  return lines.join("\n") + "\n";
}
