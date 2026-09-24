import { createContext, useCallback, useContext, useState } from "react";
import { useAppData } from "../lib/AppDataContext";
import { AboutIcon, AdmetIcon, CompareIcon, DockingIcon, DownloadIcon, LeafLattice, PredictIcon, ScreenIcon, SimilarityIcon, TargetFishingIcon, TargetInfoIcon } from "./Icons";

export type TabId =
  | "screen"
  | "predict"
  | "admet"
  | "compare"
  | "docking"
  | "target"
  | "similarity"
  | "target_fishing"
  | "downloads"
  | "about";

const TABS: { id: TabId; label: string; icon: (p: any) => JSX.Element }[] = [
  { id: "screen", label: "Screen", icon: ScreenIcon },
  { id: "predict", label: "Bioactivity Prediction", icon: PredictIcon },
  { id: "admet", label: "ADMET", icon: AdmetIcon },
  { id: "compare", label: "Compare", icon: CompareIcon },
  { id: "docking", label: "Docking", icon: DockingIcon },
  { id: "similarity", label: "Similarity", icon: SimilarityIcon },
  { id: "target_fishing", label: "Target Prediction", icon: TargetFishingIcon },
  { id: "target", label: "Target Info", icon: TargetInfoIcon },
  { id: "downloads", label: "Downloads", icon: DownloadIcon },
  { id: "about", label: "About", icon: AboutIcon },
];

export function Shell({ tab, onTab, children }: { tab: TabId; onTab: (t: TabId) => void; children: React.ReactNode }) {
  const { loading, error, targetCount } = useAppData();

  return (
    <div className="min-h-screen bg-canvas">
      <header className="sticky top-0 z-20 flex h-[58px] items-center gap-6 bg-brand-950 px-6 text-white shadow-card">
        <div className="flex items-center gap-2.5">
          <LeafLattice className="h-6 w-6 text-lime" />
          <span className="font-display text-[17px] font-semibold tracking-tight">PhytoScreen</span>
        </div>
        <nav className="flex h-full gap-1">
          {TABS.map((t) => {
            const Icon = t.icon;
            const active = tab === t.id;
            return (
              <button
                key={t.id}
                onClick={() => onTab(t.id)}
                className={`inline-flex h-full items-center gap-1.5 border-b-2 px-3.5 text-[13.5px] font-semibold transition-colors ${
                  active ? "border-lime text-white" : "border-transparent text-white/55 hover:text-white/85"
                }`}
              >
                <Icon className="h-[15px] w-[15px] opacity-90" />
                {t.label}
              </button>
            );
          })}
        </nav>
        <div className="ml-auto flex items-center gap-2 text-[12px] text-white/60">
          {error ? (
            <span className="text-clay/90">backend error</span>
          ) : loading ? (
            <span>connecting…</span>
          ) : (
            <>
              <span className="h-[7px] w-[7px] animate-pulseDot rounded-full bg-lime" />
              <b className="text-white/90">{targetCount}</b> targets loaded
            </>
          )}
        </div>
      </header>
      <main className="px-6 py-6">{children}</main>
    </div>
  );
}

export function TwoColLayout({ sidebar, children }: { sidebar: React.ReactNode; children: React.ReactNode }) {
  return (
    <div className="mx-auto grid max-w-[1800px] grid-cols-1 items-start gap-6 lg:grid-cols-[1fr_2fr]">
      <aside className="card sticky top-[78px] max-h-[calc(100vh-96px)] overflow-y-auto p-[18px]">
        {sidebar}
      </aside>
      <main className="card min-h-[60vh] overflow-hidden">{children}</main>
    </div>
  );
}

/** Sidebar open/closed is ONE shared preference across every tab, not
    per-tab state — App.tsx keeps every tab mounted at once (display:none
    toggling, so in-progress input survives switching tabs), so if each
    tab's SidebarLayout held its own useState, collapsing on "Screen"
    wouldn't visibly collapse "Docking" until a full page reload happened
    to re-read localStorage — each already-mounted instance would keep
    whatever value IT initialized to. A context provided once at the App
    root (see App.tsx's SidebarCollapseProvider) means every tab reads and
    writes the exact same state, so toggling anywhere updates everywhere
    immediately, and still remembers across reloads via localStorage. */
const SIDEBAR_STORAGE_KEY = "phytoscreen:sidebarCollapsed";
const SidebarCollapseContext = createContext<[boolean, () => void]>([false, () => {}]);

export function SidebarCollapseProvider({ children }: { children: React.ReactNode }) {
  const [collapsed, setCollapsed] = useState<boolean>(() => {
    try {
      return localStorage.getItem(SIDEBAR_STORAGE_KEY) === "1";
    } catch {
      return false;
    }
  });
  const toggle = useCallback(() => {
    setCollapsed((prev) => {
      const next = !prev;
      try {
        localStorage.setItem(SIDEBAR_STORAGE_KEY, next ? "1" : "0");
      } catch {
        /* private window / blocked storage — collapse still works for this session */
      }
      return next;
    });
  }, []);
  return <SidebarCollapseContext.Provider value={[collapsed, toggle]}>{children}</SidebarCollapseContext.Provider>;
}

export function useSidebarCollapsed(): [boolean, () => void] {
  return useContext(SidebarCollapseContext);
}

function ChevronIcon({ direction, className }: { direction: "left" | "right"; className?: string }) {
  return (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth={2.4} strokeLinecap="round" strokeLinejoin="round" className={className}>
      <path d={direction === "left" ? "M15 6l-6 6 6 6" : "M9 6l6 6-6 6"} />
    </svg>
  );
}

const TOGGLE_BTN_CLS =
  "absolute top-5 z-10 hidden h-6 w-6 items-center justify-center rounded-full border border-line bg-white text-inkmut shadow-card transition-colors hover:text-ink lg:flex";

/** The toggle button sits next to a `sticky`-positioned sidebar/content
    card. A plain `absolute` button anchored to a non-sticky ancestor does
    NOT track that card as the page scrolls — it stays at its normal-flow
    position while the card visually pins itself to the viewport, so the
    two visibly separate once there's enough scroll (barely noticeable on
    short tabs, obvious on a long one like Target Prediction's "show all"
    results). Fix: wrap the button in its own `sticky` container at the
    SAME top offset as the card (`top-[78px]`, matching every tab's own
    aside), with zero height so it adds no layout space of its own — the
    button then tracks the card exactly, since both stick at the same
    scroll-triggered point. */
function StickyToggle({ side, onClick, title }: { side: "left" | "right"; onClick: () => void; title: string }) {
  return (
    <div className="sticky top-[78px] h-0">
      <button
        type="button"
        onClick={onClick}
        title={title}
        className={`${TOGGLE_BTN_CLS} ${side === "right" ? "-right-3" : "-left-3"}`}
      >
        <ChevronIcon direction={side === "right" ? "left" : "right"} className="h-3.5 w-3.5" />
      </button>
    </div>
  );
}

/** The two-column sidebar+content grid every tab uses, with an open/close
    toggle so the sidebar can be hidden to give the content pane the full
    width — e.g. a wide results table or a big 3D view benefits from more
    room than the fixed 1fr sidebar column leaves it. `sidebar` is the
    tab's own complete <aside>...</aside> (each tab's aside already differs
    slightly — Screen/Docking use a flex column with a bottom-pinned submit
    button, others a plain scrolling card — so this wrapper only decides
    WHETHER it renders, not what's inside it) and `children` is the content
    column's own contents (rendered inside the existing "card" wrapper). */
export function SidebarLayout({
  collapsed,
  onToggle,
  sidebar,
  children,
}: {
  collapsed: boolean;
  onToggle: () => void;
  sidebar: React.ReactNode;
  children: React.ReactNode;
}) {
  // items-start (rather than the grid default, stretch) would size each
  // column to only its OWN content height -- fine when both columns are
  // similar heights, but once one column (results, e.g. Target
  // Prediction's "show all") is much taller than the other, the SHORTER
  // column's containing block ends shortly after its own short content,
  // leaving its sticky child no room to stay pinned: it unsticks and
  // scrolls away as soon as that short containing block scrolls past,
  // long before the page has scrolled through the taller column. Default
  // stretch keeps both columns at the full row height, so the sticky
  // sidebar has room to stay pinned for the whole scrollable range.
  return (
    <div className={`mx-auto grid max-w-[1800px] grid-cols-1 gap-6 ${collapsed ? "" : "lg:grid-cols-[1fr_2fr]"}`}>
      {!collapsed && (
        <div className="relative">
          <StickyToggle side="right" onClick={onToggle} title="Hide sidebar" />
          {sidebar}
        </div>
      )}
      <div className="relative min-w-0">
        {collapsed && <StickyToggle side="left" onClick={onToggle} title="Show sidebar" />}
        {children}
      </div>
    </div>
  );
}

export function SectionIntro({ title, sub }: { title: string; sub: string }) {
  return (
    <div className="mb-1">
      <h2 className="font-display text-[19px] font-medium text-ink">{title}</h2>
      <p className="mb-4 mt-0.5 text-[12.5px] text-inkmut">{sub}</p>
    </div>
  );
}

export function ResultHeader({ children }: { children: React.ReactNode }) {
  return <div className="flex flex-wrap items-center gap-x-6 gap-y-2.5 border-b border-line px-5 py-4">{children}</div>;
}
export function ResultName({ children }: { children: React.ReactNode }) {
  return <div className="font-display text-[16px] font-medium text-ink">{children}</div>;
}
export function Stat({ label, children, tip }: { label: string; children: React.ReactNode; tip?: string }) {
  return (
    <div className={`text-[12px] text-inkmut ${tip ? "cursor-help" : ""}`} title={tip}>
      {label}
      <b className="block text-[15px] font-semibold text-ink">{children}</b>
    </div>
  );
}
