"""Step 1 - pull Class A receptor list and per-residue segment annotation.

Primary source: GPCRdb. Its /residues/ endpoint returns, for every residue,
which structural segment it belongs to (TM1..TM7, ICL1-3, ECL1-3, N-term,
C-term, H8) plus the Ballesteros-Weinstein generic number. Those assignments
come from a structure-anchored profile alignment, so we do NOT have to build
our own MSA of full-length receptors (which fails badly in the loops).

Fallback: UniProt TOPO_DOM "Extracellular" features. Coarser, prediction-based,
but independent - useful as a cross-check.
"""
from __future__ import annotations

import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import quote

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from . import config as C

_SESSION = None


def session() -> requests.Session:
    global _SESSION
    if _SESSION is None:
        s = requests.Session()
        retry = Retry(total=5, backoff_factor=0.8,
                      status_forcelist=(429, 500, 502, 503, 504),
                      allowed_methods=("GET",))
        s.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=32))
        s.headers.update({"User-Agent": "gpcr-ecl-antigenic-space/1.0",
                          "Accept": "application/json"})
        _SESSION = s
    return _SESSION


def _cache_path(url: str):
    key = re.sub(r"[^A-Za-z0-9._-]", "_", url.split("//", 1)[-1])[:180]
    return C.CACHE / f"{key}.json"


def get_json(url: str, allow_404=False):
    p = _cache_path(url)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            p.unlink()
    r = session().get(url, timeout=60)
    if allow_404 and r.status_code == 404:
        p.write_text("null")
        return None
    r.raise_for_status()
    obj = r.json()
    p.write_text(json.dumps(obj))
    time.sleep(C.REQUEST_PAUSE)
    return obj


# --------------------------------------------------------------------------
# receptor list
# --------------------------------------------------------------------------
def _family_names() -> dict:
    """slug -> family name, for every descendant family of Class A."""
    fams = get_json(f"{C.GPCRDB}/proteinfamily/descendants/{C.CLASS_A_SLUG}/") or []
    return {f["slug"]: f.get("name", "") for f in fams}


def _lineage_names(slug: str, names: dict):
    parts = slug.split("_")
    out = []
    for i in range(1, len(parts) + 1):
        s = "_".join(parts[:i])
        if s in names:
            out.append(names[s])
    return out


def class_a_proteins() -> pd.DataFrame:
    fam_names = _family_names()
    prots = get_json(f"{C.GPCRDB}/proteinfamily/proteins/{C.CLASS_A_SLUG}/") or []
    rows = []
    for p in prots:
        slug = p.get("family", "") or ""
        lin = _lineage_names(slug, fam_names)
        lin_s = " | ".join(lin)
        rows.append({
            "entry_name": p.get("entry_name"),
            "accession": p.get("accession"),
            "name": p.get("name"),
            "species": p.get("species"),
            "source": p.get("source"),
            "family_slug": slug,
            "family_path": lin_s,
            "sequence": p.get("sequence") or "",
        })
    df = pd.DataFrame(rows).dropna(subset=["entry_name"])

    if C.SPECIES:
        df = df[df["species"] == C.SPECIES]
    if C.SOURCE_KEEP:
        df = df[df["source"].isin(C.SOURCE_KEEP)]

    fp = df["family_path"].fillna("")
    df["is_olfactory"] = fp.str.contains("olfact", case=False)
    df["is_taste"] = fp.str.contains("taste", case=False)
    df["is_orphan"] = fp.str.contains("orphan", case=False)
    if not C.INCLUDE_OLFACTORY:
        df = df[~df["is_olfactory"]]
    if not C.INCLUDE_TASTE:
        df = df[~df["is_taste"]]

    df = df.drop_duplicates("entry_name").sort_values("entry_name").reset_index(drop=True)
    df.to_csv(C.DATA / "proteins.csv", index=False)
    return df


# --------------------------------------------------------------------------
# per-residue segments
# --------------------------------------------------------------------------
def residues(entry_name: str):
    return get_json(f"{C.GPCRDB}/residues/{quote(entry_name)}/", allow_404=True)


def residue_table(entry_names, workers=None) -> pd.DataFrame:
    workers = workers or C.N_FETCH_WORKERS
    out, missing = [], []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for en, res in zip(entry_names, ex.map(residues, entry_names)):
            if not res:
                missing.append(en)
                continue
            for r in res:
                gn = r.get("display_generic_number") or ""
                out.append((en,
                            r.get("sequence_number"),
                            r.get("amino_acid"),
                            r.get("protein_segment"),
                            gn.split("x")[0] if gn else ""))
    df = pd.DataFrame(out, columns=["entry_name", "pos", "aa", "segment", "gn"])
    df = df.dropna(subset=["pos", "segment"])
    df["pos"] = df["pos"].astype(int)
    df.to_csv(C.DATA / "residues.csv.gz", index=False, compression="gzip")
    if missing:
        (C.DATA / "missing_residue_annotation.txt").write_text("\n".join(missing))
        print(f"[fetch] no GPCRdb residue annotation for {len(missing)} entries "
              f"(listed in data/missing_residue_annotation.txt)")
    return df


# --------------------------------------------------------------------------
# UniProt fallback / cross-check
# --------------------------------------------------------------------------
def uniprot_topology(accessions) -> pd.DataFrame:
    """Extracellular TOPO_DOM features -> N-term, ECL1, ECL2, ECL3."""
    def one(acc):
        return acc, get_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                             allow_404=True)

    rows = []
    with ThreadPoolExecutor(max_workers=C.N_FETCH_WORKERS) as ex:
        for acc, js in ex.map(one, accessions):
            if not js:
                continue
            seq = js.get("sequence", {}).get("value", "")
            ext = []
            for f in js.get("features", []):
                if f.get("type") != "Topological domain":
                    continue
                if "extracellular" not in (f.get("description") or "").lower():
                    continue
                loc = f["location"]
                s, e = loc["start"]["value"], loc["end"]["value"]
                if s and e:
                    ext.append((int(s), int(e)))
            ext.sort()
            labels = []
            for i, (s, e) in enumerate(ext):
                labels.append("N-term" if (i == 0 and s <= 3) else None)
            n_loops = 0
            for i, (s, e) in enumerate(ext):
                if labels[i] is None:
                    n_loops += 1
                    labels[i] = f"ECL{n_loops}"
            for (s, e), lab in zip(ext, labels):
                if lab in C.SEGMENTS:
                    rows.append({"accession": acc, "segment": lab,
                                 "start": s, "end": e, "seq": seq[s - 1:e]})
    df = pd.DataFrame(rows)
    df.to_csv(C.DATA / "uniprot_topology.csv", index=False)
    return df


def uniprot_annotations(accessions) -> pd.DataFrame:
    """Signal peptide end, annotated N-glycosylation positions, annotated
    disulfide pairs. GPCRdb numbering is UniProt numbering, so positions map
    directly onto the residue table."""
    def one(acc):
        return acc, get_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                             allow_404=True)

    rows = []
    with ThreadPoolExecutor(max_workers=C.N_FETCH_WORKERS) as ex:
        for acc, js in ex.map(one, accessions):
            if not js:
                continue
            sig_end, glyc, ss = 0, [], []
            for f in js.get("features", []):
                t = (f.get("type") or "").lower()
                loc = f.get("location", {})
                s = (loc.get("start") or {}).get("value")
                e = (loc.get("end") or {}).get("value")
                if t == "signal" and e:
                    sig_end = max(sig_end, int(e))
                elif t == "glycosylation" and s:
                    if "n-linked" in (f.get("description") or "").lower() \
                       or not f.get("description"):
                        glyc.append(int(s))
                elif "disulfide" in t and s and e:
                    ss.append((int(s), int(e)))
            rows.append({"accession": acc, "signal_end": sig_end,
                         "glyco_pos": ";".join(map(str, sorted(glyc))),
                         "disulfide_pairs": ";".join(f"{a}-{b}" for a, b in ss)})
    df = pd.DataFrame(rows)
    df.to_csv(C.DATA / "uniprot_annotations.csv", index=False)
    return df


def run():
    prot = class_a_proteins()
    print(f"[fetch] Class A receptors kept: {len(prot)} "
          f"({int(prot['is_orphan'].sum())} in orphan families)")
    res = residue_table(prot["entry_name"].tolist())
    print(f"[fetch] residues: {len(res):,} rows over "
          f"{res['entry_name'].nunique()} receptors")
    return prot, res


if __name__ == "__main__":
    run()
