import phronexisMark from "../assets/phronexis-mark.png";
import { LeafLattice } from "../components/Icons";

/** Content sourced from the standalone phytoscreen-about.html page (the
    project's marketing/info page), restyled entirely in this app's own
    white/green design system (canvas/surface/brand tokens, .card, the
    same font stack) instead of that page's own dark navy/teal/gold
    theme — this is an in-app tab, not a reskin of the marketing site. */

function CapabilityIcon({ path }: { path: React.ReactNode }) {
  return (
    <svg viewBox="0 0 40 40" className="h-8 w-8 stroke-brand-600" fill="none" strokeWidth={1.5}>
      {path}
    </svg>
  );
}

const CAPABILITIES = [
  {
    title: "Chemical information",
    body: "Structured compound data that grounds every natural product under study.",
    icon: (
      <>
        <circle cx="10" cy="11" r="4" />
        <circle cx="30" cy="15" r="4" />
        <circle cx="16" cy="31" r="4" />
        <line x1="10" y1="11" x2="30" y2="15" />
        <line x1="10" y1="11" x2="16" y2="31" />
        <line x1="30" y1="15" x2="16" y2="31" />
      </>
    ),
  },
  {
    title: "Molecular similarity",
    body: "Comparing structures to surface related and analogous compounds.",
    icon: (
      <>
        <circle cx="15" cy="20" r="11" />
        <circle cx="25" cy="20" r="11" />
      </>
    ),
  },
  {
    title: "Target intelligence",
    body: "Biological target context that helps ground and prioritize hypotheses.",
    icon: (
      <>
        <circle cx="20" cy="20" r="14" />
        <circle cx="20" cy="20" r="7" />
        <circle cx="20" cy="20" r="1.6" fill="currentColor" stroke="none" />
      </>
    ),
  },
  {
    title: "Computational screening",
    body: "Systematic in-silico evaluation across candidate compounds.",
    icon: (
      <>
        <rect x="6" y="6" width="28" height="28" rx="2" />
        <line x1="6" y1="16.6" x2="34" y2="16.6" />
        <line x1="6" y1="25.3" x2="34" y2="25.3" />
        <line x1="16.6" y1="6" x2="16.6" y2="34" />
        <line x1="25.3" y1="6" x2="25.3" y2="34" />
      </>
    ),
  },
];

const CONTACT_LINKS = [
  { label: "Website", href: "https://phronexis.bio", text: "phronexis.bio" },
  { label: "Email", href: "mailto:hello@phronexis.bio", text: "hello@phronexis.bio" },
  { label: "GitHub", href: "https://github.com/phronexis", text: "github.com/phronexis" },
  { label: "LinkedIn", href: "https://linkedin.com/company/phronexis", text: "linkedin.com/company/phronexis" },
];

export function AboutTab() {
  return (
    <div className="mx-auto max-w-[820px] space-y-6 pb-10">
      {/* Hero */}
      <div className="card p-8 text-center sm:p-10">
        <div className="mx-auto mb-4 flex h-14 w-14 items-center justify-center rounded-2xl bg-brand-500/10">
          <LeafLattice className="h-8 w-8 text-brand-600" />
        </div>
        <p className="text-[12.5px] font-semibold uppercase tracking-[0.08em] text-brand-700">About</p>
        <h1 className="mt-1 font-display text-[34px] font-medium text-ink sm:text-[40px]">PhytoScreen</h1>
        <p className="mt-2 font-display text-[16px] font-medium text-brand-700">
          Computational screening for natural-product discovery.
        </p>
        <p className="mx-auto mt-4 max-w-[62ch] text-[14.5px] leading-relaxed text-inkmut">
          PhytoScreen is a computational platform developed to support natural-product research and early-stage drug
          discovery. It brings together chemical information, molecular similarity, target intelligence, and
          computational screening to help researchers explore natural products and generate research hypotheses.
        </p>
      </div>

      {/* Capabilities */}
      <div className="card p-8">
        <h2 className="mb-5 font-display text-[19px] font-medium text-ink">What it brings together</h2>
        <div className="grid grid-cols-1 gap-6 sm:grid-cols-2">
          {CAPABILITIES.map((c) => (
            <div key={c.title}>
              <CapabilityIcon path={c.icon} />
              <h3 className="mt-2.5 font-display text-[15px] font-medium text-ink">{c.title}</h3>
              <p className="mt-0.5 text-[13px] text-inkmut">{c.body}</p>
            </div>
          ))}
        </div>
      </div>

      {/* Developer credit */}
      <div className="card p-8">
        <p className="mb-4 text-[11px] font-semibold uppercase tracking-[0.06em] text-inkmut">Developed by</p>
        <div className="flex flex-col items-start gap-4 sm:flex-row sm:items-center">
          <img src={phronexisMark} alt="Phronexis" className="h-11 w-auto shrink-0" />
          <div>
            <h3 className="flex flex-wrap items-center gap-2 font-display text-[17px] font-medium text-ink">
              Phronexis
              <span className="badge border border-brand-300/50 bg-brand-500/10 text-brand-700">
                Computational Biology Engineering
              </span>
            </h3>
            <p className="mt-1.5 max-w-[52ch] text-[13.5px] text-inkmut">
              Building computational tools and research infrastructure for biology, drug discovery, and life
              sciences.
            </p>
            <a
              href="https://phronexis.bio"
              target="_blank"
              rel="noopener"
              className="mt-1.5 inline-block border-b border-brand-500 text-[13px] font-medium text-brand-700 hover:text-brand-800"
            >
              phronexis.bio
            </a>
          </div>
        </div>
      </div>

      {/* Contact */}
      <div className="card p-8">
        <h2 className="font-display text-[19px] font-medium text-ink">Get in touch</h2>
        <p className="mt-1 text-[13.5px] text-inkmut">Get in touch with the Phytoscreen team.</p>
        <dl className="mt-5 divide-y divide-line">
          {CONTACT_LINKS.map((c) => (
            <div key={c.label} className="grid grid-cols-[110px_1fr] items-center py-3">
              <dt className="text-[12.5px] text-inkmut">{c.label}</dt>
              <dd>
                <a
                  href={c.href}
                  target={c.href.startsWith("http") ? "_blank" : undefined}
                  rel="noopener"
                  className="text-[13.5px] text-ink hover:text-brand-700"
                >
                  {c.text}
                </a>
              </dd>
            </div>
          ))}
        </dl>
        <a
          href="mailto:hello@phronexis.bio"
          className="mt-5 inline-flex items-center justify-center rounded-xl bg-brand-700 px-5 py-2.5 text-[14px] font-semibold text-white shadow-sm transition hover:bg-brand-800 active:scale-[0.985]"
        >
          Contact us
        </a>
      </div>

      <p className="text-center text-[12px] text-inkmut">A Phronexis project — © 2026 Phronexis. All rights reserved.</p>
    </div>
  );
}
