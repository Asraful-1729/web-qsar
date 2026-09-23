import { useEffect, useRef, useState } from "react";
import { get3Dmol } from "../lib/mol3d";

/** "Experimental vs Redocked Pose Overlay" — the visual counterpart to the
    redocking-validation RMSD number: superimposes the co-crystallized
    ligand's real experimental pose (from the PDB) against the SAME ligand
    re-docked by Vina. Only ever the reference ligand used for that
    validation, never an arbitrary submitted compound — the two poses must
    be of the identical molecule for the overlay (and the RMSD it
    illustrates) to mean anything. */
export function PoseOverlayViewer({ crystalPoseSdf, redockedPosePdb }: { crystalPoseSdf?: string | null; redockedPosePdb?: string | null }) {
  const [show, setShow] = useState(false);
  const [spin, setSpin] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const initedRef = useRef(false);

  const available = !!(crystalPoseSdf && redockedPosePdb);

  useEffect(() => {
    if (!show) {
      // The container div unmounts whenever `show` goes false (see the JSX
      // below), taking any existing 3Dmol canvas with it — resetting this
      // flag is what lets the NEXT open re-initialize into the freshly
      // remounted (empty) container, instead of the effect silently
      // no-opping because "we already initialized once" and leaving that
      // new container permanently blank (the actual bug: worked the first
      // time, stayed blank on every reopen after).
      initedRef.current = false;
      return;
    }
    if (initedRef.current || !available) return;
    initedRef.current = true;
    const $3Dmol = get3Dmol();
    if (!$3Dmol || !containerRef.current) return;
    const viewer = $3Dmol.createViewer(containerRef.current, { backgroundColor: "white" });
    viewerRef.current = viewer;
    viewer.addModel(crystalPoseSdf, "sdf");
    viewer.setStyle({ model: 0 }, { stick: { radius: 0.16, colorscheme: "grayCarbon" }, sphere: { scale: 0.22 } });
    viewer.addModel(redockedPosePdb, "pdb");
    viewer.setStyle({ model: 1 }, { stick: { radius: 0.16, colorscheme: "magentaCarbon" }, sphere: { scale: 0.22 } });
    viewer.zoomTo();
    viewer.render();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [show, available]);

  useEffect(() => {
    viewerRef.current?.spin(spin ? "y" : false);
  }, [spin]);

  if (!available) return null;

  return (
    <div className="mt-2">
      <div className="flex flex-wrap items-center gap-2">
        <button type="button" className="btn-link" onClick={() => setShow((s) => !s)}>
          Experimental vs. redocked pose overlay
        </button>
        {show && (
          <button
            type="button"
            className="btn-link"
            onClick={() => {
              viewerRef.current?.zoomTo();
              viewerRef.current?.render();
            }}
          >
            Reset view
          </button>
        )}
        {show && (
          <label className="flex cursor-pointer items-center gap-1.5 text-[12px] text-inkmut">
            <input type="checkbox" checked={spin} onChange={(e) => setSpin(e.target.checked)} />
            spin
          </label>
        )}
        {show && (
          <span className="flex items-center gap-3 text-[11px] text-inkmut">
            <span className="flex items-center gap-1">
              <span className="inline-block h-2.5 w-2.5 rounded-sm bg-[#8c8c8c]" /> experimental (crystal)
            </span>
            <span className="flex items-center gap-1">
              <span className="inline-block h-2.5 w-2.5 rounded-sm bg-[#c83cc8]" /> redocked (Vina)
            </span>
          </span>
        )}
      </div>
      {show && <div ref={containerRef} className="relative mt-2 h-[360px] w-full rounded-lg border border-line bg-white" />}
    </div>
  );
}
