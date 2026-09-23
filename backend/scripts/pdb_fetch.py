"""
Shared PDB-entry fetcher. RCSB no longer serves a legacy .pdb file for many
newer/larger depositions (mmCIF-only) — a plain files.rcsb.org/download/x.pdb
404s for those even though the structure exists. Fall back to fetching the
mmCIF and converting it with gemmi (already a project dependency), so a
receptor candidate isn't skipped just because it's a recent deposition.

That conversion has one real, unavoidable data-loss risk: the modern PDB
Chemical Component Dictionary allows 4-5 character ligand codes (e.g.
"A1AWR", confirmed on LRRK2/9C76 during a real audit — see
scripts/audit_ligand_selection.py), but legacy PDB's HETATM record has a
HARD 3-character residue-name field — gemmi's PDB writer silently
truncates ("A1AWR" -> "A1A") to fit, since there's no other way to
represent it in that 50-year-old fixed-column format. See
name_map_path()/load_name_map() below for how the true identifier is
still preserved despite that.
"""
import json
import os
import shutil
import urllib.request

import gemmi

FETCH_TIMEOUT = 30   # seconds — urlretrieve has NO timeout param of its own (only
                      # urlopen does), so an RCSB connection that stalls after accepting
                      # the socket (not a DNS/refused-connection failure, which fail fast
                      # regardless) hangs indefinitely instead of surfacing an error.


def _download(url, out_path):
    with urllib.request.urlopen(url, timeout=FETCH_TIMEOUT) as r, open(out_path, "wb") as f:
        shutil.copyfileobj(r, f)


def name_map_path(out_path):
    return out_path + ".name_map.json"


def load_name_map(out_path):
    """{truncated_3char_name: true_full_name} for any residue this
       particular fetch had to truncate to fit legacy PDB, or {} if none
       did (the common case — most components are already <=3 chars, and
       a direct legacy .pdb download never truncates anything)."""
    p = name_map_path(out_path)
    if not os.path.exists(p):
        return {}
    with open(p) as f:
        return json.load(f)


def resolve_ligand_resname(pdb_path, resname):
    """Given a ligand identifier as a caller/human understands it (e.g. the
       true "A1AWR" straight off RCSB's own ligand page), return whatever
       identifier is ACTUALLY present in pdb_path's HETATM records —
       Bio.PDB reads resname directly off those fixed columns, so a 4-5
       character CCD code that fetch_pdb()'s gemmi fallback had to
       truncate to fit them (see module docstring) won't text-match
       against the file unless translated back to that truncated form
       first. A no-op (returns resname unchanged) whenever resname is
       already <=3 characters, or no truncation happened for this fetch,
       or resname doesn't match anything in the sidecar map — those are
       exactly the cases where no translation is needed or none is
       possible, not errors."""
    if resname is None or len(resname) <= 3:
        return resname
    name_map = load_name_map(pdb_path)
    for truncated, full in name_map.items():
        if full == resname:
            return truncated
    return resname


def fetch_pdb(pdb_id, out_path):
    """Writes a legacy-format PDB file to out_path, converting from mmCIF if
       RCSB has no legacy .pdb for this entry. If that conversion truncates
       any residue name, writes a sidecar name_map_path(out_path) recording
       {truncated: true_full_name} — load it with load_name_map() to
       recover what a truncated HETATM code in THIS file actually is;
       receptor_prep.py's resolve_ligand_resname() uses it so callers can
       keep passing the real, human-meaningful identifier without needing
       to know this quirk exists."""
    try:
        _download(f"https://files.rcsb.org/download/{pdb_id}.pdb", out_path)
        return out_path
    except urllib.error.HTTPError as e:
        if e.code != 404:
            raise
    cif_path = out_path + ".cif"
    _download(f"https://files.rcsb.org/download/{pdb_id}.cif", cif_path)
    structure = gemmi.read_structure(cif_path)
    structure.setup_entities()

    name_map = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if len(res.name) > 3:
                    name_map[res.name[:3]] = res.name
        break   # first model only — matches receptor_prep.py's own convention

    structure.write_pdb(out_path)
    os.remove(cif_path)
    if name_map:
        with open(name_map_path(out_path), "w") as f:
            json.dump(name_map, f)
    return out_path
