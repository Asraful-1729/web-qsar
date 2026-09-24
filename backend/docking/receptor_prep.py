"""
Receptor preparation (build-time, once per target).  [WRITTEN â VALIDATE ON YOUR MACHINE]

raw PDB (with co-crystallised ligand)
  -> extract reference ligand            (Bio.PDB)
  -> strip waters/HETATMs                (Bio.PDB)
  -> repair: add missing atoms + H       (PDBFixer / OpenMM)
  -> grid box from reference ligand      (grid_box_from_ligand â TESTED math)
  -> receptor PDBQT                       (Meeko mk_prepare_receptor)
  -> save DockingProfile to docking_registry.json

!! I could not execute this end-to-end (PDBFixer/real PDBs unavailable
   in development). The pure-Python parts (ligand extraction, box math, .gpf text)
   are tested; the external-tool steps are standard but UNVALIDATED. VALIDATE by
   re-docking the co-crystallised ligand and confirming RMSD < ~2 A before trusting
   any result for this target.

Dependencies (install on your machine):
   conda install -c conda-forge pdbfixer openmm
   pip install meeko                      # provides mk_prepare_receptor
"""
import os, json, shutil, subprocess, tempfile
import numpy as np

from .profile import grid_box_from_ligand, REGISTRY, registry_lock, write_registry_json
from subprocess_util import hidden_subprocess_kwargs

# The additive/cryoprotectant/cofactor blacklist and MW floor used to decide
# "is this HETATM group a real ligand" live in scripts/select_receptor.py
# (imported lazily below, matching this module's existing convention) —
# this used to have its own smaller, drifted COMMON_ADDITIVES copy (19
# entries vs. select_receptor's 90+, missing e.g. HEM/NAD/FAD/ATP), which
# meant extract_reference_ligand() and list_ligands() could each decide
# something different was "the ligand" for the exact same structure. One
# shared blacklist now, everywhere a HETATM group needs judging.

# Sanity-check floors/ceilings for the ligand extract_reference_ligand()
# picks — real inhibitors bound in a real pocket fall well inside these;
# anything outside gets FLAGGED (see ligand_sanity_warnings in
# build_receptor()), never silently rejected or silently trusted.
MIN_HEAVY_ATOMS = 5        # matches list_ligands()'s own floor for "not a stray fragment/ion"
MAX_LIGAND_MW = 2000.0     # generous — this app screens natural products, some legitimately large
MIN_POCKET_RESIDUES = 3    # fewer than this means the "ligand" barely touches the protein at all

# Average atomic weights (Da), heavy atoms only — real depositions almost
# never have hydrogens, so this is deliberately a heavy-atom-only proxy
# (~85-95% of true MW for typical drug-like ligands), not a precise value.
# Good enough to catch "this is obviously way too small/large," which is
# all it's used for; anything unrecognised falls back to a carbon-ish
# average rather than raising.
_ATOMIC_WEIGHTS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "P": 30.974, "S": 32.06, "CL": 35.45, "BR": 79.904, "I": 126.904,
    "B": 10.811, "SI": 28.085, "SE": 78.971, "NA": 22.99, "K": 39.098,
    "MG": 24.305, "CA": 40.078, "ZN": 65.38, "FE": 55.845, "MN": 54.938,
}


def _approx_heavy_atom_mw(atoms):
    """Sum of average atomic weights over an iterable of Bio.PDB Atoms —
       see _ATOMIC_WEIGHTS above for what this is/isn't good for."""
    total = 0.0
    for atom in atoms:
        el = (getattr(atom, "element", "") or "").strip().upper()
        if not el:
            nm = atom.get_name().strip()
            el = nm[0].upper() if nm else ""
        total += _ATOMIC_WEIGHTS.get(el, 12.0)
    return total


def _count_pdb_atoms(path):
    """Cheap line-scan atom count (ATOM+HETATM) — used only to report real
       before/after facts in build_receptor()'s prep_report, not for
       anything that needs Bio.PDB's actual structure model."""
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                n += 1
    return n


def _count_pdbqt_atoms(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                n += 1
    return n


# ---------- step 1: parse + extract reference ligand ----------
def extract_reference_ligand(pdb_path, ref_resname=None, chain=None, resnum=None):
    """Return (ref_ligand_coords Nx3, ref_resname, ref_atoms, ref_mw) for the
       co-crystallised inhibitor. If ref_resname is None, pick the largest
       non-additive HETATM group — see MIN_HEAVY_ATOMS/MAX_LIGAND_MW above
       for the sanity floors build_receptor() checks the result against
       (this function itself never rejects a pick on size grounds; when
       ref_resname is explicit, that's a deliberate human/curated choice
       that should surface as a warning if it looks odd, not get silently
       overridden here).
       resnum: only needed to disambiguate the rare case of two separate
       copies of the SAME ligand resname within the SAME chain (distinct
       chains already disambiguate on their own) — see list_ligands(),
       which is what a caller picking a specific candidate from a
       multi-ligand structure uses to get an exact (resname, chain,
       resnum) triple in the first place.

       ref_resname may be given as the true CCD identifier (e.g. "A1AWR")
       even if pdb_fetch.py's gemmi fallback had to truncate it to fit
       legacy PDB's 3-character HETATM field to actually get written into
       pdb_path — see pdb_fetch.resolve_ligand_resname. The returned
       ref_resname is always the truncated, FILE-CONSISTENT form (what a
       later text search against this same file, e.g. locate_ligand_near,
       must use), never the true one; build_receptor() separately records
       the caller's original (possibly-untruncated) argument for display."""
    from Bio.PDB import PDBParser
    from scripts.select_receptor import ADDITIVE_BLACKLIST
    from scripts.pdb_fetch import resolve_ligand_resname
    ref_resname = resolve_ligand_resname(pdb_path, ref_resname)
    s = PDBParser(QUIET=True).get_structure("x", pdb_path)
    candidates = {}
    for model in s:
        for ch in model:
            if chain and ch.id != chain:
                continue
            for res in ch:
                het = res.id[0].strip()               # '' for standard, 'H_XXX'/'W' for hetero
                if not het:
                    continue
                name = res.resname.strip()
                if name in ADDITIVE_BLACKLIST:
                    continue
                if ref_resname and name != ref_resname:
                    continue
                if resnum is not None and res.id[1] != resnum:
                    continue
                atoms = list(res.get_atoms())
                coords = np.array([a.coord for a in atoms], float)
                key = (name, ch.id, res.id[1])
                candidates[key] = (coords, atoms)
        break                                          # first model only
    if not candidates:
        raise ValueError("no reference ligand found (check ref_resname / additives list)")
    # largest heavy-atom group wins (deterministic)
    key = max(candidates, key=lambda k: len(candidates[k][0]))
    coords, atoms = candidates[key]
    mw = _approx_heavy_atom_mw(atoms)
    return coords, key[0], coords.shape[0], mw


def list_ligands(pdb_path, min_heavy_atoms=5, max_peptide_gap=2, min_peptide_run=3):
    """Every real (non-additive, non-covalent-modification) small-molecule
       ligand co-crystallized in this raw structure — a real crystal
       structure often has MULTIPLE distinct bound ligands (different
       copies of the same inhibitor in separate chains, a genuinely
       different molecule in a second site, or an unrelated fragment), and
       extract_reference_ligand() above always silently picks just the
       single LARGEST one as THE reference. This lists all of them so a
       caller can offer 'redock centered on a different one' instead of
       being stuck with that one automatic pick.

       Three filters, each catching a different kind of non-ligand HETATM
       group that a naive 'every hetero residue' scan would wrongly offer
       as a selectable ligand:
       1. select_receptor.py's ADDITIVE_BLACKLIST (buffers, ions,
          cryoprotectants, glycosylation sugars, nucleotide-analog
          cofactors) — the same one extract_reference_ligand() above uses.
       2. Bio.PDB's is_aa(..., standard=False): catches a single modified
          amino acid embedded mid-polypeptide-chain (e.g. SEP/TPO/CSO —
          phosphorylated/oxidized residues, common on kinase activation
          loops) — a real part of the PROTEIN, not a ligand to dock
          against, even though PDB format marks it HETATM.
       3. Contiguous-run grouping: some structures instead have a whole
          PEPTIDE/macrocycle inhibitor built from non-standard (often
          D-)amino acids Bio.PDB doesn't recognize via #2 (e.g. a p53-MDM2
          D-peptide inhibitor, DAL/DCY/DGL/... one HETATM per residue,
          resnum 1-17 in sequence) — is_aa() alone won't catch these, but
          >=`min_peptide_run` hetero residues in the SAME chain with
          resnums within `max_peptide_gap` of each other is a reliable
          structural signature of 'this is one polymer chain, not several
          independent small molecules,' regardless of what it's built
          from. This pipeline docks single small molecules via Vina from a
          SMILES string, so a whole peptide inhibitor is out of scope for
          'pick an alternate ligand to redock against' either way — the
          WHOLE run is excluded rather than offering its 17 residues as if
          they were 17 different candidate pockets.

       A small heavy-atom-count floor (default 5) additionally drops
       leftover fragments/monatomic ions the blacklist doesn't name.
       Distinct (resname, chain, resnum) copies of the SAME real ligand
       (e.g. one per protein chain in a crystallographic dimer) are listed
       separately — they sit in different pockets, so which one becomes
       the reference genuinely matters for where the docking box lands."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa
    from scripts.select_receptor import ADDITIVE_BLACKLIST
    s = PDBParser(QUIET=True).get_structure("x", pdb_path)

    by_chain = {}
    for model in s:
        for ch in model:
            for res in ch:
                het = res.id[0].strip()
                if not het:
                    continue
                name = res.resname.strip()
                if name in ADDITIVE_BLACKLIST:
                    continue
                if is_aa(name, standard=False):
                    continue
                coords = np.array([a.coord for a in res.get_atoms()], float)
                if coords.shape[0] < min_heavy_atoms:
                    continue
                centroid = coords.mean(axis=0)
                by_chain.setdefault(ch.id, []).append({
                    "resname": name, "chain": ch.id, "resnum": res.id[1], "n_atoms": int(coords.shape[0]),
                    "center": [round(float(v), 3) for v in centroid],
                })
        break   # first model only, same convention as extract_reference_ligand

    out = []
    for chain_id, residues in by_chain.items():
        residues.sort(key=lambda r: r["resnum"])
        run = [residues[0]]
        for r in residues[1:]:
            if r["resnum"] - run[-1]["resnum"] <= max_peptide_gap:
                run.append(r)
            else:
                if len(run) < min_peptide_run:
                    out.extend(run)
                run = [r]
        if len(run) < min_peptide_run:
            out.extend(run)

    out.sort(key=lambda r: -r["n_atoms"])
    return out


# ---------- step 2: strip to protein only ----------
def strip_to_protein(pdb_path, out_pdb, chain=None):
    """Standard residues only (no HETATM/water). If chain is given, keep ONLY
       that chain — many PDB entries deposit 2+ copies of the same protein in
       the asymmetric unit (e.g. a crystallographic dimer), and merging them
       into one 'receptor' both docks against a physically wrong target (two
       overlapping copies) and can confuse bond-perception in atoms close to
       the chain-chain interface (RDKit inferring spurious cross-chain bonds)."""
    from Bio.PDB import PDBParser, PDBIO, Select

    class ProteinOnly(Select):
        def accept_residue(self, res):
            if res.id[0] != " ":
                return False
            if chain and res.get_parent().id != chain:
                return False
            return True

    s = PDBParser(QUIET=True).get_structure("x", pdb_path)
    io = PDBIO(); io.set_structure(s)
    io.save(out_pdb, ProteinOnly())
    return out_pdb


# ---------- step 3: repair (PDBFixer) ----------
def _minimize_added_atoms(fixer, max_iterations=300):
    """PDBFixer's addMissingAtoms()/addMissingHydrogens() place new atoms
       from geometric templates with no clash-checking against the rest of
       the residue — for a residue whose crystal structure was missing a
       sidechain (or part of one), this can drop the rebuilt atoms almost on
       top of one another. Observed on ADRB1/7BU6: GLN284 and ARG366 each got
       a sidechain CG landing 1.56 A from their own NE2/CZ — a 1-3
       (two-bonds-apart) distance that should be ~2.3-2.4 A. Meeko's
       receptor-to-mol step infers bonds from interatomic distance alone
       (see _drop_oxt below for the other flavor of this bug), so that clash
       reads as a real bond and blows the atom's valence, crashing
       mk_prepare_receptor.py with an RDKit AtomValenceException.

       A short vacuum energy minimization relaxes exactly these local
       clashes — real bonded geometry already sits in a deep energy well, so
       nothing conformationally meaningful moves; confirmed on the case
       above (energy dropped from +4e5 to -5e4 kJ/mol, CG-NE2 distance
       corrected to 2.45 A, ~2s on the CPU platform for a ~7200-atom
       receptor) before the structure ever reaches Meeko."""
    from openmm import app, unit, LocalEnergyMinimizer, VerletIntegrator, Context, Platform
    forcefield = app.ForceField("amber14-all.xml")
    system = forcefield.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff,
                                     constraints=None, rigidWater=False)
    integrator = VerletIntegrator(1.0 * unit.femtoseconds)
    try:
        platform = Platform.getPlatformByName("CPU")
    except Exception:
        platform = Platform.getPlatformByName("Reference")
    context = Context(system, integrator, platform)
    context.setPositions(fixer.positions)
    LocalEnergyMinimizer.minimize(context, maxIterations=max_iterations)
    fixer.positions = context.getState(getPositions=True).getPositions()


def repair_receptor(pdb_in, pdb_out, ph=7.0):
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except Exception as e:
        raise RuntimeError(f"PDBFixer/OpenMM not installed: {e}. "
                           f"conda install -c conda-forge pdbfixer openmm")
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    try:
        _minimize_added_atoms(fixer)
    except Exception:
        pass   # best-effort clash relief — a forcefield-template mismatch here must not block receptor prep
    with open(pdb_out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    _drop_oxt(pdb_out)
    return pdb_out


def _drop_oxt(pdb_path):
    """Remove C-terminal OXT atoms in place. PDBFixer completes every chain
       terminus with one (standard chemistry), but Meeko's PDB->mol step
       perceives bond orders from interatomic distance alone: the terminal
       carbon's two near-equidistant C-O contacts (O and OXT) both get read
       as double bonds, pushing its valence to 5 and crashing receptor prep
       on EVERY chain terminus. Vina scores a rigid receptor purely from atom
       positions/types (no formal bond orders), so dropping this one terminal
       oxygen costs nothing chemically relevant to docking — confirmed fix
       against CHEMBL1862_ABL1/1IEP (chains A and B both hit this)."""
    with open(pdb_path) as f:
        lines = f.readlines()
    kept = [ln for ln in lines
            if not (ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "OXT")]
    with open(pdb_path, "w") as f:
        f.writelines(kept)


def _mk_prepare_receptor_inprocess(clean_pdb, stem):
    """Call meeko's mk_prepare_receptor CLI script's own main() in-process,
       instead of shelling out to the 'mk_prepare_receptor.py' console-
       script launcher pip generates in a venv's Scripts/ folder — that
       launcher genuinely doesn't exist in a frozen PyInstaller build (no
       venv, no PATH-discoverable executable), even though the underlying
       meeko PACKAGE is fully bundled (--collect-all meeko in
       BUILD_WINDOWS.md — the CLI script's own source, meeko/cli/
       mk_prepare_receptor.py, is a normal file inside that package).
       Sidesteps subprocess+PATH entirely — works identically whether
       frozen or not, and surfaces real Python exceptions (an RDKit
       sanitization failure, say) directly instead of parsed subprocess
       stderr text. Returns True if meeko ran (success or a real prep
       failure, raised as RuntimeError); False if meeko itself can't be
       imported at all, so the caller can fall back to prepare_receptor."""
    import sys, io, contextlib
    try:
        from meeko.cli import mk_prepare_receptor as _mk_cli
    except Exception:
        return False
    old_argv = sys.argv
    sys.argv = ["mk_prepare_receptor.py", "--read_pdb", clean_pdb, "-o", stem, "-p"]
    err_buf = io.StringIO()
    try:
        with contextlib.redirect_stderr(err_buf):
            _mk_cli.main()
    except SystemExit as e:
        if e.code not in (0, None):
            raise RuntimeError(f"mk_prepare_receptor failed on {clean_pdb}:\n{err_buf.getvalue()}") from e
    except Exception as e:
        raise RuntimeError(f"mk_prepare_receptor failed on {clean_pdb}: {e}\n{err_buf.getvalue()}") from e
    finally:
        sys.argv = old_argv
    return True


# ---------- step 4: receptor PDBQT (Meeko) ----------
def receptor_to_pdbqt(clean_pdb, out_pdbqt):
    """Prefer Meeko's mk_prepare_receptor (in-process — see
       _mk_prepare_receptor_inprocess); fall back to prepare_receptor
       (ADFR/MGLTools, a real external binary the desktop app doesn't
       bundle — genuinely optional)."""
    stem = out_pdbqt[:-6] if out_pdbqt.endswith(".pdbqt") else out_pdbqt
    if _mk_prepare_receptor_inprocess(clean_pdb, stem):
        cand = stem + ".pdbqt"
        if os.path.exists(cand):
            if cand != out_pdbqt:
                shutil.move(cand, out_pdbqt)
            return out_pdbqt
    if shutil.which("prepare_receptor"):
        subprocess.run(["prepare_receptor", "-r", clean_pdb, "-o", out_pdbqt],
                       check=True, capture_output=True, text=True, timeout=120, **hidden_subprocess_kwargs())
        return out_pdbqt
    raise RuntimeError("no receptor-prep tool found. `pip install meeko` (mk_prepare_receptor) "
                       "or install ADFR/MGLTools prepare_receptor.")


# ---------- binding site (pocket residues, for showing WHY the box sits
# where it does — not just an opaque center/size) ----------
def largest_protein_chain(pdb_path):
    """The chain id with the most standard-residue (polymer) atoms — a
       best-effort single-chain pick for build_receptor()'s no-ligand
       fallback below, where there's no co-crystallized ligand to infer
       'which chain actually matters' from (see strip_to_protein's own
       docstring for why restricting to ONE chain matters: many real
       depositions carry 2+ copies of the same protein). Returns None if
       the structure has no standard-residue atoms at all (never expected
       in practice — strip_to_protein's own downstream failure would be
       the more informative error at that point)."""
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("x", pdb_path)
    counts = {}
    for model in s:
        for ch in model:
            n = sum(1 for res in ch if res.id[0] == " ")
            if n:
                counts[ch.id] = n
        break
    if not counts:
        return None
    return max(counts, key=counts.get)


def locate_ligand_near(pdb_path, resname, near_point, chain=None):
    """Among every HETATM group named `resname` in pdb_path (any chain unless
       `chain` is given), return the heavy-atom coords of whichever copy's
       centroid is closest to `near_point`. Chain isn't persisted in
       docking_registry.json, but the box center a profile actually produced
       IS — distance-matching against that recovers the exact ligand copy a
       profile was built from even for structures depositing 2+ copies in
       the asymmetric unit (e.g. 1IEP/ABL1, see detect_chain.py)."""
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("x", pdb_path)
    near = np.asarray(near_point, float)
    best = None
    for model in s:
        for ch in model:
            if chain and ch.id != chain:
                continue
            for res in ch:
                if res.resname.strip() != resname:
                    continue
                coords = np.array([a.coord for a in res.get_atoms()], float)
                d = float(np.linalg.norm(coords.mean(axis=0) - near))
                if best is None or d < best[0]:
                    best = (d, coords)
        break
    if best is None:
        raise ValueError(f"ligand '{resname}' not found in {pdb_path}")
    return best[1]


def pocket_residues(clean_pdb_path, ref_coords, cutoff=5.0):
    """Receptor residues with any atom within `cutoff` A of any reference-
       ligand atom — the binding site shown to users before docking runs,
       instead of leaving 'why does the box sit here' implicit."""
    from Bio.PDB import PDBParser, NeighborSearch
    s = PDBParser(QUIET=True).get_structure("x", clean_pdb_path)
    ns = NeighborSearch(list(s.get_atoms()))
    seen = {}
    for pt in ref_coords:
        for atom in ns.search(pt, cutoff):
            res = atom.get_parent()
            key = (res.get_parent().id, res.id[1], res.resname.strip())
            seen[key] = True
    out = [{"chain": k[0], "resnum": k[1], "resname": k[2]} for k in seen]
    out.sort(key=lambda r: (r["chain"], r["resnum"]))
    return out


def all_residues(clean_pdb_path):
    """Every standard (non-heteroatom) residue in the receptor — the full
       amino acid sequence, structured per-residue (chain/resnum/resname)
       rather than a flat one-letter string, since that's what a manual
       binding-site picker actually needs: something to list, highlight in
       3D, and select. Unlike pocket_residues() this isn't centered on any
       reference ligand at all, so it works even for a target with no
       co-crystallized ligand and needs no ref_coords — 'manual' binding-
       site definition means picking from the WHOLE protein, not just
       refining the already automatically-detected pocket neighborhood."""
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("x", clean_pdb_path)
    out = []
    for model in s:
        for ch in model:
            for res in ch:
                if res.id[0] != " ":   # skip HETATM/water — same filter strip_to_protein() uses
                    continue
                out.append({"chain": ch.id, "resnum": res.id[1], "resname": res.resname.strip()})
        break
    out.sort(key=lambda r: (r["chain"], r["resnum"]))
    return out


def box_from_residues(clean_pdb_path, residues, padding=8.0, min_size=20.0):
    """center/box_size (grid_box_from_ligand's contract) from a user-picked
       residue subset — Advanced Settings' 'define a custom binding site
       from residues' path, an alternative to typing raw coordinates."""
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("x", clean_pdb_path)
    want = {(r["chain"], int(r["resnum"])) for r in residues}
    coords = []
    for model in s:
        for ch in model:
            for res in ch:
                if (ch.id, res.id[1]) in want:
                    coords.extend(a.coord for a in res.get_atoms())
        break
    if not coords:
        raise ValueError("none of the given residues were found in this receptor")
    return grid_box_from_ligand(coords, padding=padding, min_size=min_size)


def box_from_receptor(clean_pdb_path, padding=4.0, min_size=20.0):
    """Whole-protein bounding box — the 'blind docking' path (no pocket
       assumption at all, as opposed to the ligand-centered or residue-
       selected site-specific boxes above). Smaller padding than the
       site-specific default (4 A vs 8 A): the box already spans the whole
       receptor, so there's no need for the extra margin a small, focused
       pocket box wants. Vina still has to search a MUCH larger volume than
       a site-specific box, so this is inherently slower and less reliable
       per-site than a validated pocket — callers should surface that as a
       caveat, not silently treat blind results as equally trustworthy."""
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure("x", clean_pdb_path)
    coords = [a.coord for a in s.get_atoms()]
    if not coords:
        raise ValueError(f"no atoms found in {clean_pdb_path}")
    return grid_box_from_ligand(coords, padding=padding, min_size=min_size)


# ---------- orchestration ----------
def build_receptor(pdb_path, target_id, name=None, ref_resname=None, chain=None, resnum=None,
                   out_dir="docking_targets", padding=8.0, progress=None, ligand_chain=None):
    """Runs the full strip -> repair -> PDBQT pipeline and returns a profile
       dict with ABSOLUTE file paths — no docking_registry.json I/O. Used by
       both prepare_receptor() (build-time, persists into the shared
       registry below) and app.py's on-demand 'Advanced Settings' manual
       structure override (deliberately never persisted — an expert's
       per-request pick shouldn't silently overwrite the vetted default, and
       staying out of the shared registry file avoids racing a concurrent
       batch_validate.py run that owns writes to it).

       resnum: passed straight through to extract_reference_ligand() — see
       its docstring; only needed to disambiguate two same-resname copies
       within the same chain, which list_ligands() callers can supply.

       ligand_chain: which chain the REFERENCE LIGAND's own HETATM records
       are filed under, if different from `chain` (which protein chain the
       RECEPTOR gets stripped to). These are usually the same chain, but
       not always: in a hetero-complex deposition, a ligand sitting at a
       chain-chain interface is filed under whichever chain letter the
       depositor happened to assign it, which doesn't have to be the chain
       that actually forms most of its binding contacts — confirmed on
       GENE_UBA2/6XOG (ligand VAY filed under chain C, the complex's SUMO1
       subunit, but 86 of its 93 nearby protein atoms are on chain B, the
       real UBA2/SAE2 subunit — chain B is what should be built as "the
       UBA2 receptor") and GENE_PDK3/1Y8P (RED filed under chain B, a
       different protein entirely — DLAT's E2 component — while 39 of 43
       contacts are on chain A, the real PDK3 protein). Defaults to `chain`
       when omitted, which is correct for the ordinary case (a ligand
       filed under the same chain it actually binds).

       progress, if given, is called with a short human-readable label
       before each real stage starts — the on-demand manual-structure path
       (app.py's /api/docking/receptor/custom) surfaces these live so a
       user watching a ~1-2 minute prep isn't just staring at one static
       'preparing...' message the whole time."""
    if ligand_chain is None:
        ligand_chain = chain
    from scripts.select_receptor import MIN_LIGAND_MW

    def _p(label):
        if progress:
            progress(label)

    tdir = os.path.join(out_dir, target_id)
    os.makedirs(tdir, exist_ok=True)
    # B13 — real, non-fabricated before/after facts for each stage, surfaced
    # to the caller as `prep_report` below: every number here comes from an
    # actual line-count of the file that stage just produced, not a label
    # restating what the stage does (which is all `progress()`'s callback
    # gives the UI live — this is the same pipeline's DURABLE record of it).
    prep_report = []
    n_input_atoms = _count_pdb_atoms(pdb_path)

    _p("Extracting reference ligand")
    try:
        ref_coords, ref_name, n_ref, ref_mw = extract_reference_ligand(pdb_path, ref_resname, ligand_chain, resnum=resnum)
        center, box_size = grid_box_from_ligand(ref_coords, padding=padding)
        prep_report.append({"label": "Extract reference ligand", "detail": f"{ref_name}: {n_ref} atom(s), ~{ref_mw:.0f} Da"})
    except ValueError:
        if ref_resname:
            raise   # an explicit pick (curated data or a user's own choice) that doesn't exist here is a real error — surface it, don't silently fall back
        # ref_resname was never given (automatic pick) and NOTHING in this
        # structure qualifies as a real small-molecule ligand — e.g.
        # PLEC/2ODV, whose only HETATM group is a cryoprotectant (see
        # scripts/audit_ligand_selection.py's no_reference_ligand_recorded
        # targets; ~70 more genes across panel_results_v2.csv hit this same
        # wall). That doesn't make the STRUCTURE unusable — it just means
        # there's no pocket to center a site-specific box on, so this
        # builds a Blind-only receptor instead of failing outright.
        ref_coords = ref_name = None
        n_ref, ref_mw = 0, 0.0
        center = box_size = None
        if chain is None:
            chain = largest_protein_chain(pdb_path)
        prep_report.append({"label": "Extract reference ligand", "detail":
                            "no real small-molecule ligand found (only buffers/cryoprotectants/ions, if anything) — "
                            "building a Blind-only receptor, no site-specific default"})

    _p("Stripping to protein-only (removing waters/heteroatoms)")
    prot = strip_to_protein(pdb_path, os.path.join(tdir, "protein_raw.pdb"), chain=chain)
    n_stripped_atoms = _count_pdb_atoms(prot)
    prep_report.append({"label": "Strip to protein-only", "detail":
                        f"{n_input_atoms} → {n_stripped_atoms} atoms "
                        f"({n_input_atoms - n_stripped_atoms} water(s)/heteroatom(s) removed)"})

    _p("Repairing (PDBFixer: missing atoms/residues, hydrogenation)")
    clean = repair_receptor(prot, os.path.join(tdir, "receptor_clean.pdb"))
    n_repaired_atoms = _count_pdb_atoms(clean)
    prep_report.append({"label": "Repair (PDBFixer)", "detail":
                        f"{n_stripped_atoms} → {n_repaired_atoms} atoms "
                        f"({n_repaired_atoms - n_stripped_atoms} missing atom(s)/hydrogen(s) added)"})

    _p("Building PDBQT (Meeko: atom typing, charges)")
    rec_pdbqt = receptor_to_pdbqt(clean, os.path.join(tdir, "receptor.pdbqt"))
    n_pdbqt_atoms = _count_pdbqt_atoms(rec_pdbqt)
    prep_report.append({"label": "Build PDBQT (Meeko)", "detail":
                        f"{n_pdbqt_atoms} atom(s) typed with AutoDock atom types + partial charges "
                        "(receptor kept rigid for docking — no torsions here)"})

    _p("Computing binding site (pocket residues, grid box)")
    try:
        binding_site_residues = pocket_residues(clean, ref_coords, cutoff=5.0) if ref_coords is not None else []
    except Exception:
        binding_site_residues = []   # non-fatal — box/docking still work without this display data
    try:
        all_residues_list = all_residues(clean)
    except Exception:
        all_residues_list = []   # non-fatal — manual site-picking just won't have a full list to show
    try:
        blind_center, blind_box_size = box_from_receptor(clean)
    except Exception:
        blind_center, blind_box_size = None, None   # non-fatal — blind mode just won't be offered for this structure
    if center and box_size:
        prep_report.append({"label": "Compute binding site", "detail":
                            f"{len(binding_site_residues)} pocket residue(s) within 5.0 Å of the reference ligand "
                            f"({len(all_residues_list)} total residue(s) in the receptor); "
                            f"box center ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}), "
                            f"size {box_size[0]:.1f} × {box_size[1]:.1f} × {box_size[2]:.1f} Å"})
    else:
        prep_report.append({"label": "Compute binding site", "detail":
                            f"no site-specific box ({len(all_residues_list)} total residue(s) in the receptor); "
                            + (f"blind (whole-protein) box size {blind_box_size[0]:.1f} × {blind_box_size[1]:.1f} × {blind_box_size[2]:.1f} Å"
                               if blind_box_size else "blind box unavailable too")})

    # Flag, never silently reject or silently trust: the additive blacklist
    # already rules out named non-ligands, but it can't be exhaustive (a
    # cofactor/fragment it doesn't name yet, or an explicit ref_resname
    # that turns out to be a poor choice) — these are the residual checks
    # a human should look at before trusting this target's default box.
    # None of this applies when there's no ligand at all (the automatic-
    # pick-found-nothing case above) — that's an already-understood,
    # expected state (site_source below says so directly), not something
    # to flag for review the way an unusually small/large PICKED ligand is.
    ligand_sanity_warnings = []
    if ref_coords is not None:
        if n_ref < MIN_HEAVY_ATOMS:
            ligand_sanity_warnings.append(
                f"only {n_ref} heavy atom(s) in '{ref_name}' — unusually small for a real inhibitor "
                f"(floor: {MIN_HEAVY_ATOMS})")
        if ref_mw < MIN_LIGAND_MW:
            ligand_sanity_warnings.append(
                f"'{ref_name}' is ~{ref_mw:.0f} Da, below the {MIN_LIGAND_MW:.0f} Da floor — "
                "may be an additive/fragment the blacklist doesn't name")
        if ref_mw > MAX_LIGAND_MW:
            ligand_sanity_warnings.append(
                f"'{ref_name}' is ~{ref_mw:.0f} Da, unusually large — check this isn't a mis-flagged "
                "peptide/chain rather than a single ligand")
        if len(binding_site_residues) < MIN_POCKET_RESIDUES:
            ligand_sanity_warnings.append(
                f"only {len(binding_site_residues)} receptor residue(s) within 5 Å of '{ref_name}' — "
                f"it may not actually sit in a real binding pocket (floor: {MIN_POCKET_RESIDUES})")
        if ligand_sanity_warnings:
            prep_report.append({"label": "Ligand sanity check", "detail": "; ".join(ligand_sanity_warnings)})

    return {
        "prep_report": prep_report,
        "all_residues": all_residues_list,
        "target_id": target_id, "name": name or target_id,
        "pdb_source": os.path.basename(pdb_path), "reference_ligand_resname": ref_name,
        # Set only when the caller's own ref_resname argument (the TRUE CCD
        # id, e.g. "A1AWR") had to be resolved to a different, truncated
        # in-file form (ref_name, e.g. "A1A") — see extract_reference_ligand
        # / pdb_fetch.resolve_ligand_resname. reference_ligand_resname stays
        # the truncated form deliberately (locate_ligand_near does a live
        # text search against the raw file by that exact string), so this
        # is purely an informational/display field recovering what a
        # truncated 3-char code in THIS registry entry actually means.
        "reference_ligand_full_id": ref_resname if (ref_resname and ref_resname != ref_name) else None,
        "chain": chain,                # persisted so a later revert/repair (see batch_validate.py's
                                        # _accept_or_revert) can rebuild an IDENTICAL receptor from
                                        # scratch — without this, a revert only restored the registry's
                                        # center/box_size/pdb_source, not which chain was stripped to,
                                        # so a rebuild from raw pdb_source alone could silently include
                                        # extra chains never present in the originally-validated receptor.
        # Only set when it differs from `chain` — the ligand's own HETATM
        # chain label, when that's a DIFFERENT chain than the one the
        # receptor was stripped to (see the ligand_chain param docstring
        # above). Needed for the same revert/rebuild-from-scratch reason as
        # `chain` itself: without it, a rebuild would default ligand_chain
        # back to `chain` and fail to find the ligand at all.
        "ligand_chain": ligand_chain if ligand_chain != chain else None,
        "receptor_pdbqt": os.path.abspath(rec_pdbqt),
        "receptor_pdb": os.path.abspath(clean),
        "center": center, "box_size": box_size,
        "binding_site_residues": binding_site_residues,
        "blind_center": blind_center, "blind_box_size": blind_box_size,
        "site_source": "co-crystal_ligand" if ref_name else "none_validated" if ref_coords is None else "auto",
        "reference_ligand_mw": round(ref_mw, 1) if ref_coords is not None else None,
        "ligand_sanity_warnings": ligand_sanity_warnings,
        # "crystallographic_annotation": the caller told us exactly which
        # HETATM group is the real ligand (from panel_results_v2.csv, a
        # curated pick, or a user's own Advanced Settings choice).
        # "automatic": nobody knew, so this fell back to "largest
        # non-additive HETATM group" — buildable ("ok" doesn't fail) is not
        # the same claim as biologically validated, which is exactly what
        # build_status distinguishes below.
        # "none": no ref_resname given AND nothing qualified — see the
        # no-ligand fallback above; this target is Blind-only, same as the
        # audit's no_reference_ligand_recorded targets, just discovered
        # on-demand here instead of during the original batch build.
        "ligand_selection_method": "crystallographic_annotation" if ref_resname else ("none" if ref_coords is None else "automatic"),
        "build_status": "review_required" if ligand_sanity_warnings else "ok",
    }


def prepare_receptor(pdb_path, target_id, name=None, ref_resname=None, chain=None,
                     out_dir="docking_targets", padding=8.0):
    """Build-time onboarding: build_receptor() + persist into the shared
       docking_registry.json with portable (basename) paths."""
    profile = build_receptor(pdb_path, target_id, name=name, ref_resname=ref_resname,
                             chain=chain, out_dir=out_dir, padding=padding)
    # store PORTABLE basenames, not absolute paths — resolved at load time
    # relative to DOCKING_TARGETS_DIR/<target_id>/ (see profile.load_profile),
    # so the registry keeps working after the project moves or is packaged.
    profile = dict(profile, receptor_pdbqt=os.path.basename(profile["receptor_pdbqt"]),
                   receptor_pdb=os.path.basename(profile["receptor_pdb"]))
    with registry_lock():   # see registry_lock's docstring — needed once multiple targets can validate concurrently
        reg = {}
        if os.path.exists(REGISTRY):
            data = json.load(open(REGISTRY))
            reg = {t["target_id"]: t for t in data.get("targets", [])}
        reg[target_id] = profile
        write_registry_json({"targets": list(reg.values())})
    return profile