import { useCallback, useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { Modal } from "./Modal";
import type { AdvancedDockingState } from "../lib/useAdvancedDocking";
import { residueKey } from "../lib/useAdvancedDocking";
import {
  applyProteinStyle,
  clearSurfaces,
  drawBoxShapes,
  fetchTextCached,
  get3Dmol,
  HANDLE_SPECS,
  MIN_BOX_DIM,
  type ProteinStyle,
} from "../lib/mol3d";
import { StyleSelect } from "./StyleSelect";

type HandleKind = "move" | "resize";
interface HandleVM {
  key: string;
  kind: HandleKind;
  axis?: number;
  sign?: number;
  color?: string;
  x: number;
  y: number;
  visible: boolean;
}

export function BindingSiteModal({
  adv,
  targetId,
  onClose,
}: {
  adv: AdvancedDockingState;
  targetId: string;
  onClose: () => void;
}) {
  const { site, dockingMode, selected, toggleResidue, boxOverride, setBoxOverrideFromDrag, effectiveBox } = adv;
  const blind = dockingMode === "blind";
  const [residueFilter, setResidueFilter] = useState("");

  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const [ready, setReady] = useState(false);
  const [gizmoOn, setGizmoOn] = useState(false);
  const [handles, setHandles] = useState<HandleVM[]>([]);
  const [readout, setReadout] = useState<{ center: number[]; size: number[] } | null>(null);
  const [style, setStyle] = useState<ProteinStyle>("cartoon");
  const surfaceIdsRef = useRef<number[]>([]);

  /** The receptor's base representation, plus whatever residue highlight
      is currently selected layered on top — shared by the full rebuild,
      the cheap re-highlight, and the style switcher so all three agree
      on what "the current look" is instead of drifting out of sync. */
  const applyReceptorLook = useCallback(
    (viewer: any, s: ProteinStyle) => {
      if (!site?.receptorUrl) return;
      surfaceIdsRef.current = clearSurfaces(viewer, surfaceIdsRef.current);
      applyProteinStyle(viewer, { model: 0 }, s, (ids) => (surfaceIdsRef.current = ids));
      if (s === "cartoon") viewer.setStyle({ model: 0 }, { cartoon: { color: "lightgrey", opacity: 0.5 } });
      const pool = site.allResidues?.length ? site.allResidues : site.residues;
      const selResi = [...new Set(pool.filter((r) => selected.has(residueKey(r))).map((r) => r.resnum))];
      if (selResi.length) viewer.addStyle({ model: 0, resi: selResi }, { stick: { radius: 0.18, colorscheme: "yellowCarbon" } });
    },
    [site, selected]
  );

  const positionHandles = useCallback(() => {
    const viewer = viewerRef.current;
    const [center, size] = effectiveBox();
    if (!viewer || !gizmoOn || blind || !center || !size) {
      setHandles([]);
      return;
    }
    const specs: { key: string; kind: HandleKind; axis?: number; sign?: number; color?: string }[] = [
      ...HANDLE_SPECS.map((h) => ({ key: h.key, kind: "resize" as const, axis: h.axis, sign: h.sign, color: h.color })),
      { key: "move", kind: "move" },
    ];
    // Fires on every 3Dmol view-change too (setViewChangeCallback below),
    // including mid-rotation when the user clicks-drags on the model
    // itself rather than precisely on a handle — modelToScreen can return
    // a transient null for a frame or two while the camera matrix is
    // mid-update. Falling back to {x:0,y:0,visible:false} there made
    // handles visibly vanish/flicker during that rotation, which is the
    // "can not hold it, it disappears" complaint — enabling drag editing
    // should keep the handles up until the user turns it back off, not
    // whenever the camera merely moves. Keeping each handle's last known
    // good screen position instead means a transient failed projection is
    // invisible to the user rather than a flicker.
    setHandles((prev) => {
      const prevByKey = new Map(prev.map((h) => [h.key, h]));
      return specs.map((h) => {
        let pt: { x: number; y: number; z: number };
        if (h.kind === "move") pt = { x: center[0], y: center[1], z: center[2] };
        else {
          pt = { x: center[0], y: center[1], z: center[2] };
          const axisKey = (["x", "y", "z"] as const)[h.axis!];
          pt[axisKey] += h.sign! * (size[h.axis!] / 2);
        }
        const s = viewer.modelToScreen(pt);
        if (!s) {
          const last = prevByKey.get(h.key);
          return last ? { ...h, x: last.x, y: last.y, visible: last.visible } : { ...h, x: 0, y: 0, visible: false };
        }
        return { ...h, x: s.x, y: s.y, visible: true };
      });
    });
  }, [effectiveBox, gizmoOn, blind]);

  const updateReadout = useCallback(() => {
    const [center, size] = effectiveBox();
    if (!center || !size) {
      setReadout(null);
      return;
    }
    setReadout({ center, size });
  }, [effectiveBox]);

  // positionHandles is a fresh function every time gizmoOn/effectiveBox/blind
  // change, but viewer.setViewChangeCallback below is only REGISTERED once
  // per receptor load (the full-rebuild effect only re-runs on a new site) —
  // that registration closes over whatever positionHandles WAS at mount
  // time, permanently, with gizmoOn=false baked in (its value when the
  // modal first opens, before the user has checked "Enable drag editing").
  // Every later camera change (rotating/zooming the model, not just
  // dragging a handle) invoked that stale closure, which always took the
  // `!gizmoOn` early-exit and wiped every handle via setHandles([]) — this
  // was the actual "handles disappear the moment I touch the model, can't
  // hold them" bug, not a transient projection failure. Routing the
  // registered callback through a ref that's kept current means it always
  // calls whichever positionHandles closure is live right now.
  const positionHandlesRef = useRef(positionHandles);
  useEffect(() => {
    positionHandlesRef.current = positionHandles;
  }, [positionHandles]);

  // ---- full rebuild: new site / target / structure ----
  useEffect(() => {
    let cancelled = false;
    setReady(false);
    async function run() {
      const $3Dmol = get3Dmol();
      if (!$3Dmol || !containerRef.current || !site) return;
      containerRef.current.innerHTML = "";
      const viewer = $3Dmol.createViewer(containerRef.current, { backgroundColor: "white" });
      viewerRef.current = viewer;
      viewer.setViewChangeCallback(() => positionHandlesRef.current());

      let modelIdx = 0;
      let resiList: number[] = [];
      surfaceIdsRef.current = [];
      if (site.receptorUrl) {
        const pdb = await fetchTextCached(site.receptorUrl);
        if (cancelled) return;
        if (pdb) {
          viewer.addModel(pdb, "pdb");
          applyReceptorLook(viewer, style);
          resiList = [...new Set(site.residues.map((r) => r.resnum))];
          modelIdx++;
        }
      }
      if (site.ligandUrl) {
        try {
          const sdf = await fetchTextCached(site.ligandUrl);
          if (cancelled) return;
          if (sdf) {
            viewer.addModel(sdf, "sdf");
            viewer.setStyle({ model: modelIdx }, { stick: { radius: 0.18, colorscheme: "cyanCarbon" } });
            modelIdx++;
          }
        } catch {
          /* optional */
        }
      }
      const [center, size] = effectiveBox();
      if (center && size) drawBoxShapes(viewer, center, size);
      // Start focused on the automatically-detected pocket (if any) —
      // every residue is still listed and selectable, so the user can
      // freely pan/zoom out to pick from elsewhere (e.g. an allosteric
      // site); this just avoids opening on a zoomed-out, unfocused view
      // for the common case where the automatic pocket is exactly right.
      if (resiList.length) viewer.zoomTo({ model: 0, resi: resiList });
      else viewer.zoomTo();
      viewer.render();
      if (cancelled) return;
      setReady(true);
      positionHandles();
      updateReadout();
    }
    run();
    return () => {
      cancelled = true;
      viewerRef.current = null;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [site?.receptorUrl, site?.ligandUrl]);

  // ---- cheap residue re-highlight (no rebuild) ----
  useEffect(() => {
    const viewer = viewerRef.current;
    // `ready` gates this: the full-rebuild effect sets viewerRef.current
    // synchronously (before its own `await fetchTextCached(...)`), so a
    // viewer can exist here with zero models still added — model 0 not
    // existing yet crashes 3Dmol's setStyle (TypeError deep inside
    // GLViewer.applyToModels), which — uncaught in an effect — unmounts the
    // whole React tree, not just this modal. Only touch models once the
    // rebuild has actually finished adding them.
    if (!viewer || !site || !ready) return;
    applyReceptorLook(viewer, style);
    viewer.render();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selected, ready, style]);

  // ---- box-only redraw: mode / override changes ----
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !ready) return;
    const [center, size] = effectiveBox();
    viewer.removeAllShapes();
    if (center && size) drawBoxShapes(viewer, center, size);
    viewer.render();
    positionHandles();
    updateReadout();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [dockingMode, boxOverride, ready]);

  useEffect(() => {
    positionHandles();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [gizmoOn]);

  useEffect(() => {
    const onResize = () => positionHandles();
    window.addEventListener("resize", onResize);
    window.addEventListener("scroll", onResize, true);
    return () => {
      window.removeEventListener("resize", onResize);
      window.removeEventListener("scroll", onResize, true);
    };
  }, [positionHandles]);

  const dragState = useRef<{ handle: HandleVM; startX: number; startY: number; center0: number[]; size0: number[]; raf: number | null } | null>(null);

  const applyDrag = useCallback(
    (dx: number, dy: number) => {
      const d = dragState.current;
      const viewer = viewerRef.current;
      if (!d || !viewer) return;
      const delta = viewer.screenOffsetToModel(dx, dy);
      const dv = [delta.x, delta.y, delta.z];
      let center = d.center0.slice();
      let size = d.size0.slice();
      if (d.handle.kind === "move") {
        center = [center[0] + dv[0], center[1] + dv[1], center[2] + dv[2]];
      } else {
        const a = d.handle.axis!;
        const s = d.handle.sign!;
        let newSize = d.size0[a] + s * dv[a];
        newSize = Math.max(MIN_BOX_DIM, newSize);
        const actualDelta = newSize - d.size0[a];
        center[a] = d.center0[a] + (s * actualDelta) / 2;
        size[a] = newSize;
      }
      const c3: [number, number, number] = [center[0], center[1], center[2]];
      const s3: [number, number, number] = [size[0], size[1], size[2]];
      setBoxOverrideFromDrag(c3, s3);
      viewer.removeAllShapes();
      drawBoxShapes(viewer, c3, s3);
      viewer.render();
      positionHandles();
      setReadout({ center: c3, size: s3 });
    },
    [setBoxOverrideFromDrag, positionHandles]
  );

  const startDrag = useCallback(
    (handle: HandleVM, e: React.MouseEvent) => {
      e.preventDefault();
      e.stopPropagation();
      const [center0, size0] = effectiveBox();
      if (!center0 || !size0) return;
      dragState.current = { handle, startX: e.clientX, startY: e.clientY, center0, size0, raf: null };
      const onMove = (mev: MouseEvent) => {
        const d = dragState.current;
        if (!d) return;
        if (d.raf) return;
        d.raf = requestAnimationFrame(() => {
          if (dragState.current) dragState.current.raf = null;
          applyDrag(mev.clientX - d.startX, mev.clientY - d.startY);
        });
      };
      const onUp = () => {
        document.removeEventListener("mousemove", onMove);
        document.removeEventListener("mouseup", onUp);
        dragState.current = null;
      };
      document.addEventListener("mousemove", onMove);
      document.addEventListener("mouseup", onUp);
    },
    [effectiveBox, applyDrag]
  );

  // Always the WHOLE receptor's residues — every amino acid is listed, not
  // just the automatically-detected pocket, so nothing is hidden from the
  // user; the pocket residues just arrive pre-checked (see useAdvancedDocking's
  // loadBindingSite), and any residue can be added or removed from there.
  const residues = site?.allResidues?.length ? site.allResidues : site?.residues || [];
  const selectedCount = residues.filter((r) => selected.has(residueKey(r))).length;
  const filteredResidues = residueFilter.trim()
    ? residues.filter((r) => {
        const q = residueFilter.trim().toLowerCase();
        return (
          r.resname.toLowerCase().includes(q) ||
          String(r.resnum).includes(q) ||
          r.chain.toLowerCase().includes(q) ||
          `${r.chain}:${r.resnum}`.toLowerCase().includes(q)
        );
      })
    : residues;

  return (
    <Modal
      title={`Binding site — ${targetId}${blind ? " (blind — whole protein)" : ""}`}
      sub={
        blind
          ? "Whole-protein search box — fixed, not editable. Switch to Site-specific mode to move/resize it."
          : "Drag the orange handle to move the box, the colored handles to resize it."
      }
      onClose={onClose}
    >
      <div className="flex h-full">
        <div className="relative flex-1 bg-white" ref={containerRef}>
          <div className="absolute right-1.5 top-1.5 z-10">
            <StyleSelect value={style} onChange={setStyle} />
          </div>
        </div>
        <div className="w-[240px] shrink-0 overflow-y-auto border-l border-line bg-canvas/70 p-3.5">
          {!blind && (
            <label className="mb-3 flex cursor-pointer items-center gap-2 text-[13px] font-medium text-ink">
              <input type="checkbox" checked={gizmoOn} onChange={(e) => setGizmoOn(e.target.checked)} />
              Enable drag editing
            </label>
          )}
          <p className="field-hint mb-3">
            {blind
              ? "Blind mode searches the entire receptor surface — the box is fixed to the whole-protein bounding volume and not user-editable."
              : "Orange = move box. Red/green/blue = resize X/Y/Z (anchored on the opposite face)."}
          </p>
          <label className="field-label">Current box</label>
          <div className="rounded-lg border border-line bg-surface px-2.5 py-2 text-[12px] text-ink">
            {readout ? (
              <>
                <div className="flex justify-between py-px">
                  <span className="text-inkmut">Center (Å)</span>
                  <b>{readout.center.map((v) => v.toFixed(1)).join(", ")}</b>
                </div>
                <div className="flex justify-between py-px">
                  <span className="text-inkmut">Size (Å)</span>
                  <b>{readout.size.map((v) => v.toFixed(1)).join(" × ")}</b>
                </div>
              </>
            ) : (
              <span className="text-inkmut">—</span>
            )}
          </div>
          <label className="field-label">
            {`All residues (${residues.length}) ${
              blind ? "— informational only in blind mode" : "— pocket residues are pre-checked; toggle any to add/remove them"
            }`}
          </label>
          {residues.length > 15 && (
            <input
              className="field-input mb-1.5 text-[12.5px]"
              placeholder="Filter by residue, number, or chain (e.g. LEU, 145, A:145)"
              value={residueFilter}
              onChange={(e) => setResidueFilter(e.target.value)}
            />
          )}
          <div className="max-h-[220px] overflow-y-auto rounded-lg border border-line">
            {!residues.length && <div className="p-2 text-[12.5px] text-inkmut">No residue data for this structure.</div>}
            {!!residues.length && !filteredResidues.length && (
              <div className="p-2 text-[12.5px] text-inkmut">No residues match "{residueFilter}".</div>
            )}
            {filteredResidues.map((r) => {
              const key = residueKey(r);
              return (
                <label key={key} className="flex cursor-pointer items-center gap-2 border-b border-line/70 px-2.5 py-1.5 text-[12.5px] last:border-0 hover:bg-surface2/60">
                  <input type="checkbox" checked={selected.has(key)} onChange={(e) => toggleResidue(key, e.target.checked)} />
                  <span className={selected.has(key) ? "font-semibold text-amber" : "text-inkmut"}>
                    {r.resname} {r.resnum} · chain {r.chain}
                  </span>
                </label>
              );
            })}
          </div>
          {!!residues.length && !blind && (
            <div className="field-hint">
              {selectedCount === 0
                ? "No residues selected — pick at least one to build a binding box."
                : `Binding box set from ${selectedCount} of ${residues.length} residue(s).`}
            </div>
          )}
        </div>
      </div>
      {!blind &&
        gizmoOn &&
        handles
          .filter((h) => h.visible)
          .map((h) =>
            createPortal(
              <div
                key={h.key}
                onMouseDown={(e) => startDrag(h, e)}
                className="absolute z-[9999] rounded-full border-2 border-white shadow-md"
                style={{
                  left: h.x,
                  top: h.y,
                  width: h.kind === "move" ? 20 : 16,
                  height: h.kind === "move" ? 20 : 16,
                  marginLeft: h.kind === "move" ? -10 : -8,
                  marginTop: h.kind === "move" ? -10 : -8,
                  background: h.kind === "move" ? "#f59e0b" : h.color,
                  borderRadius: h.kind === "move" ? 6 : 999,
                  cursor: h.kind === "move" ? "move" : "grab",
                }}
              />,
              document.body
            )
          )}
    </Modal>
  );
}
