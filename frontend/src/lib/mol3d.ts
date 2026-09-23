// Small shared helpers for the two 3Dmol.js integrations (binding-site
// viewer w/ drag gizmo, and the docking-result ball-and-stick pose viewer).
// 3Dmol is loaded globally via the vendored <script> tag (see index.html) —
// there's no npm package dependency, matching the original app.

declare global {
  interface Window {
    $3Dmol: any;
  }
}

export function get3Dmol(): any | null {
  return typeof window !== "undefined" && window.$3Dmol ? window.$3Dmol : null;
}

const _textCache = new Map<string, Promise<string | null>>();
export function fetchTextCached(url: string): Promise<string | null> {
  if (!_textCache.has(url)) {
    _textCache.set(
      url,
      fetch(url)
        .then((r) => (r.ok ? r.text() : null))
        .catch(() => null)
    );
  }
  return _textCache.get(url)!;
}

/** Wireframe-only outline — deliberately no filled/translucent box. A
    filled box at any opacity > 0 is still a solid mesh WebGL has to
    depth-sort against everything else in the scene (the receptor sticks,
    the ligand); 3Dmol draws shapes in insertion order rather than true
    back-to-front sorting, so atoms that happen to fall "behind" the box's
    triangles from the current camera angle silently fail their depth test
    and don't render — until the camera moves and the sort order changes,
    which is exactly the "some atoms inside the box don't show until I
    rotate" symptom this fixes. A pure wireframe has no fill triangles to
    occlude anything, so nothing to depth-sort against in the first place. */
export function drawBoxShapes(viewer: any, center: [number, number, number], size: [number, number, number]) {
  if (!center || !size) return;
  viewer.addBox({
    center: { x: center[0], y: center[1], z: center[2] },
    dimensions: { w: size[0], h: size[1], d: size[2] },
    color: "red",
    wireframe: true,
    linewidth: 2,
  });
}

/** Combines a receptor PDB and a docked ligand-pose PDB into one
    downloadable "complex" file — same two structures PoseViewer already
    loads as separate 3Dmol models for VISUALIZATION, just concatenated
    as text for a real file. Strips any existing END/ENDMDL from each
    part (so the receptor's own terminator doesn't cut the file short)
    and adds a TER between chains plus a single trailing END. */
export function combinePdbText(receptorPdb: string, posePdb: string): string {
  const strip = (s: string) =>
    s
      .split("\n")
      .filter((line) => !/^(END|ENDMDL)\s*$/.test(line.trim()))
      .join("\n")
      .replace(/\n+$/, "");
  return `${strip(receptorPdb)}\nTER\n${strip(posePdb)}\nEND\n`;
}

export const MIN_BOX_DIM = 4.0;

/** Protein representation styles offered by the style switcher in every
    3D viewer (ReceptorPreview, PoseViewer, BindingSiteModal). "surface"
    variants and "mesh" go through addSurface — heavier than the atom/
    bond styles, and need their own cleanup (see clearSurfaces) since
    3Dmol tracks them separately from setStyle. */
export type ProteinStyle = "cartoon" | "stick" | "line" | "sphere" | "surface" | "surfaceHydrophobicity" | "mesh";

export const PROTEIN_STYLE_OPTIONS: { value: ProteinStyle; label: string }[] = [
  { value: "cartoon", label: "Ribbon (cartoon)" },
  { value: "stick", label: "Stick" },
  { value: "line", label: "Line" },
  { value: "sphere", label: "Sphere" },
  { value: "surface", label: "Surface" },
  { value: "surfaceHydrophobicity", label: "Hydrophobicity surface" },
  { value: "mesh", label: "Mesh" },
];

// Kyte & Doolittle hydrophobicity scale (most positive = most hydrophobic).
const KD_SCALE: Record<string, number> = {
  ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5, MET: 1.9, ALA: 1.8,
  GLY: -0.4, THR: -0.7, SER: -0.8, TRP: -0.9, TYR: -1.3, PRO: -1.6,
  HIS: -3.2, GLU: -3.5, GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5,
};

function lerp(a: number, b: number, t: number): number {
  return a + (b - a) * t;
}

/** Blue (hydrophilic) -> white -> orange (hydrophobic), same convention
    PyMOL/Chimera use for this coloring. */
function hydrophobicityColorForValue(v: number): string {
  const t = (v + 4.5) / 9; // -4.5..4.5 -> 0..1
  const lo = { r: 0x2b, g: 0x6c, b: 0xb0 }; // blue
  const mid = { r: 0xff, g: 0xff, b: 0xff }; // white
  const hi = { r: 0xe2, g: 0x7d, b: 0x2a }; // orange
  const [a, b, u] = t < 0.5 ? [lo, mid, t / 0.5] : [mid, hi, (t - 0.5) / 0.5];
  const r = Math.round(lerp(a.r, b.r, u));
  const g = Math.round(lerp(a.g, b.g, u));
  const bl = Math.round(lerp(a.b, b.b, u));
  return `#${[r, g, bl].map((x) => x.toString(16).padStart(2, "0")).join("")}`;
}

// A pre-sampled color ramp for 3Dmol's addSurface `map: {prop, gradient}`
// coloring path (see applyProteinStyle's surfaceHydrophobicity case) —
// deliberately NOT a `colorfunc` closure. 3Dmol computes surface geometry
// across several Web Workers; a colorfunc can't cross that boundary (JS
// closures aren't structured-cloneable), so it runs on the MAIN THREAD once
// per vertex as each worker posts geometry back — thousands of synchronous
// calls for a full receptor, and worse, clearSurfaces() only deletes the
// rendered geometry, not the still-running workers, so switching styles
// while one is mid-flight just queues more blocking work behind the click
// instead of cancelling it. `map`/`gradient` is plain, clonable data (a
// property name + numeric min/max + a color array), so the worker resolves
// vertex colors entirely on its own — nothing runs on the main thread per
// vertex, and there's nothing left to cancel.
//
// 3Dmol's CustomLinear gradient buckets N colors into N equal-width
// segments interpolating colors[i]->colors[i+1], with the LAST bucket
// always flat (no i+1 to interpolate toward) — not a naive N-1-segment
// gradient. Oversampling the smooth blue-white-orange curve at many stops
// keeps that flat final bucket negligibly small instead of visibly
// clipping the top of the range to solid orange.
const HYDROPHOBICITY_GRADIENT_STOPS = 20;
const HYDROPHOBICITY_GRADIENT_COLORS: string[] = Array.from({ length: HYDROPHOBICITY_GRADIENT_STOPS }, (_, i) =>
  hydrophobicityColorForValue(-4.5 + (9 * i) / (HYDROPHOBICITY_GRADIENT_STOPS - 1))
);

/** Removes every surface this viewer is currently tracking — call before
    applying a non-surface style, and before adding a new surface (3Dmol
    layers surfaces rather than replacing them, so switching from one
    surface style to another without this leaves the old one behind). */
export function clearSurfaces(viewer: any, ids: number[]): number[] {
  for (const id of ids) {
    try {
      viewer.removeSurface(id);
    } catch {
      /* already gone */
    }
  }
  return [];
}

/** Applies one of PROTEIN_STYLE_OPTIONS to `sel` (a 3Dmol selection spec,
    e.g. {model: 0}). Non-surface styles resolve synchronously; surface
    styles are async in 3Dmol (addSurface returns before the mesh is
    built) — pass the returned surface id(s) to clearSurfaces() later via
    onSurfaceIds, and call viewer.render() again once it resolves (3Dmol
    does this internally too, but an explicit render right after keeps
    the underlying atom style visible immediately instead of blank). */
export function applyProteinStyle(viewer: any, sel: any, style: ProteinStyle, onSurfaceIds?: (ids: number[]) => void): void {
  switch (style) {
    case "cartoon":
      viewer.setStyle(sel, { cartoon: { color: "lightgrey" } });
      break;
    case "stick":
      viewer.setStyle(sel, { stick: { radius: 0.15, colorscheme: "grayCarbon" } });
      break;
    case "line":
      viewer.setStyle(sel, { line: {} });
      break;
    case "sphere":
      viewer.setStyle(sel, { sphere: { scale: 0.3, colorscheme: "grayCarbon" } });
      break;
    case "surface": {
      viewer.setStyle(sel, { cartoon: { color: "lightgrey" } });
      const id = viewer.addSurface("VDW", { opacity: 0.85, color: "white" }, sel);
      onSurfaceIds?.([id]);
      break;
    }
    case "surfaceHydrophobicity": {
      viewer.setStyle(sel, { cartoon: { color: "lightgrey" } });
      // Stamp each atom's KD hydrophobicity as a plain numeric property
      // (0 = neutral, for ligand atoms/waters/anything unrecognised) so
      // addSurface's `map` coloring can read it worker-side — see the
      // HYDROPHOBICITY_GRADIENT_COLORS comment above for why this replaces
      // a colorfunc closure instead of just wrapping one.
      for (const atom of viewer.selectedAtoms(sel)) {
        if (!atom.properties) atom.properties = {};
        atom.properties.hydrophob = KD_SCALE[(atom.resn || "").toUpperCase()] ?? 0;
      }
      const id = viewer.addSurface(
        "VDW",
        {
          opacity: 0.9,
          map: { prop: "hydrophob", gradient: { gradient: "linear", min: -4.5, max: 4.5, colors: HYDROPHOBICITY_GRADIENT_COLORS } },
        },
        sel
      );
      onSurfaceIds?.([id]);
      break;
    }
    case "mesh": {
      viewer.setStyle(sel, { cartoon: { color: "lightgrey" } });
      const id = viewer.addSurface("SES", { opacity: 1, wireframe: true, color: "#6b7280" }, sel);
      onSurfaceIds?.([id]);
      break;
    }
  }
}

/** Shows every HETATM group (ligands, ions, cofactors, crystallographic
    waters) as ball-and-stick, layered on top of whatever protein style is
    active. Needed because none of PROTEIN_STYLE_OPTIONS render heteroatoms
    on their own — 3Dmol's "cartoon" (and every other polymer-backbone
    style) only draws residues that are part of a recognized protein/
    nucleic-acid chain, so a HETATM group gets assigned that style and
    simply renders as nothing. Used for the "before" (raw, as-deposited)
    receptor preview so it actually shows what RCSB shows — the original
    file's real ligands/waters — not just a bare protein ribbon; the
    "after" (stripped/repaired) side has none of this left to show. */
export function applyHeteroStyle(viewer: any, sel: any): void {
  // Waters (single, unbonded oxygens) fall back to the sphere and render
  // as their default CPK red — plenty visible on their own. Everything
  // else (the actual ligand(s), ions, cofactors) gets a magenta carbon
  // scheme specifically because it reads clearly against both the white
  // background and the light-grey cartoon, unlike a neutral/white one.
  viewer.setStyle({ ...sel, hetflag: true }, { stick: { radius: 0.18, colorscheme: "magentaCarbon" }, sphere: { scale: 0.25 } });
}

export const HANDLE_SPECS = [
  { key: "x1", axis: 0, sign: 1, color: "#ef4444" },
  { key: "x-1", axis: 0, sign: -1, color: "#ef4444" },
  { key: "y1", axis: 1, sign: 1, color: "#22c55e" },
  { key: "y-1", axis: 1, sign: -1, color: "#22c55e" },
  { key: "z1", axis: 2, sign: 1, color: "#3b82f6" },
  { key: "z-1", axis: 2, sign: -1, color: "#3b82f6" },
] as const;
