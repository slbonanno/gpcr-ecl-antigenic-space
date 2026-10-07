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

import numpy as np
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
        retry = Retry(total=6, backoff_factor=1.5,
                      status_forcelist=(429, 500, 502, 503, 504),
                      allowed_methods=("GET",),
                      respect_retry_after_header=True,
                      raise_on_status=False)
        s.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=32))
        s.headers.update({"User-Agent": "gpcr-ecl-antigenic-space/1.0",
                          "Accept": "application/json"})
        _SESSION = s
    return _SESSION


def _cache_path(url: str):
    key = re.sub(r"[^A-Za-z0-9._-]", "_", url.split("//", 1)[-1])[:180]
    return C.CACHE / f"{key}.json"


def get_json(url: str, allow_404=False, tolerant=False, pause=None):
    """Cached GET.

    tolerant=True returns None instead of raising when the server is simply
    unavailable (503/429 after retries, connection reset, read timeout). A
    transient failure is NEVER written to the cache, so re-running picks it up.
    A genuine 404 with allow_404 IS cached, because it will stay a 404.
    """
    p = _cache_path(url)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            p.unlink()
    try:
        r = session().get(url, timeout=60)
    except requests.RequestException as e:
        if tolerant:
            return None
        raise
    if allow_404 and r.status_code == 404:
        p.write_text("null")
        return None
    if r.status_code >= 400:
        if tolerant:
            return None
        r.raise_for_status()
    try:
        obj = r.json()
    except ValueError:
        if tolerant:
            return None
        raise
    p.write_text(json.dumps(obj))
    time.sleep(C.REQUEST_PAUSE if pause is None else pause)
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
    nm = df["name"].fillna("")
    df["is_olfactory"] = fp.str.contains("olfact", case=False)
    df["is_taste"] = fp.str.contains("taste", case=False)

    # Orphan status. GPCRdb does not label every deorphanisation-pending
    # receptor with the word "orphan" in its family name, so a family-path
    # match alone silently returns zero. Fall back on the naming convention:
    # UniProt calls a receptor with no assigned ligand "G-protein coupled
    # receptor <n>" / "Probable G-protein coupled receptor <n>", and GPCRdb
    # entry names for those are gpr<n>_human / gprc5a_human etc.
    by_family = fp.str.contains("orphan", case=False)
    by_name = nm.str.contains(r"^(?:probable |putative )?g[- ]?protein[- ]coupled receptor",
                              case=False, regex=True)
    by_entry = df["entry_name"].str.match(r"^gpr(?:c)?\d", case=False, na=False)
    df["is_orphan"] = by_family | by_name | by_entry
    df["orphan_evidence"] = np.select(
        [by_family, by_name | by_entry], ["family", "naming"], default="")

    if by_family.sum() == 0:
        print("[fetch] note: no GPCRdb family is named 'orphan'; orphan flag "
              "fell back to receptor naming convention. Family names are "
              "dumped to data/family_tree.csv — check them if the count "
              "looks wrong.")
    pd.DataFrame(sorted(fam_names.items()), columns=["slug", "name"]).to_csv(
        C.DATA / "family_tree.csv", index=False)
    if not C.INCLUDE_OLFACTORY:
        df = df[~df["is_olfactory"]]
    if not C.INCLUDE_TASTE:
        df = df[~df["is_taste"]]

    df = df.drop_duplicates("entry_name").sort_values("entry_name").reset_index(drop=True)
    df.to_csv(C.DATA / "proteins.csv", index=False)

    top = (df["family_path"].fillna("").str.split(" | ", regex=False)
             .str[1].fillna("(unclassified)").value_counts())
    print("[fetch] receptors by top-level family:")
    print(top.to_string())
    return df


# --------------------------------------------------------------------------
# per-residue segments
# --------------------------------------------------------------------------
def residues(entry_name: str):
    return get_json(f"{C.GPCRDB}/residues/{quote(entry_name)}/",
                    allow_404=True, tolerant=True)


def residue_table(entry_names, workers=None) -> pd.DataFrame:
    workers = workers or C.N_FETCH_WORKERS

    def safe(en):
        try:
            return residues(en)
        except Exception:
            return None

    out, missing = [], []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for en, res in zip(entry_names, ex.map(safe, entry_names)):
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


def _uniprot_batches(accessions, size=None):
    """UniProt's per-accession endpoint will 503 you off the service if you hit
    it a few hundred times in parallel. Use the stream endpoint instead and ask
    for many accessions at once: ~3 requests instead of ~300."""
    size = size or C.UNIPROT_BATCH
    accs = [a for a in accessions if isinstance(a, str) and a]
    for i in range(0, len(accs), size):
        chunk = accs[i:i + size]
        q = quote(" OR ".join(f"accession:{a}" for a in chunk))
        url = ("https://rest.uniprot.org/uniprotkb/stream"
               f"?query={q}&format=json"
               "&fields=accession,ft_signal,ft_carbohyd,ft_disulfid")
        yield chunk, url


def uniprot_annotations(accessions) -> pd.DataFrame:
    """Signal peptide end, annotated N-glycosylation positions, annotated
    disulfide pairs. GPCRdb numbering is UniProt numbering, so positions map
    straight onto the residue table.

    Entirely optional: if UniProt is unreachable this returns an empty frame
    and the pipeline carries on without signal-peptide trimming.
    """
    entries, failed = [], 0
    for chunk, url in _uniprot_batches(accessions):
        js = get_json(url, tolerant=True, pause=C.UNIPROT_PAUSE)
        if js is None:
            failed += len(chunk)
            continue
        entries.extend(js.get("results", js) if isinstance(js, dict) else js)

    if failed:
        print(f"[fetch] UniProt unreachable for {failed} accessions; "
              "re-run later to fill them in (nothing was cached)")
    if not entries:
        print("[fetch] no UniProt annotations - signal peptides will NOT be "
              "trimmed and sequons will be used in place of annotated glycans")
        return pd.DataFrame(columns=["accession", "signal_end", "glyco_pos",
                                     "disulfide_pairs"])

    rows = []
    for js in entries:
        acc = js.get("primaryAccession") or js.get("accession")
        if not acc:
            continue
        sig_end, glyc, ss = 0, [], []
        for f in js.get("features", []):
            ftype = (f.get("type") or "").lower()
            loc = f.get("location", {})
            s = (loc.get("start") or {}).get("value")
            e = (loc.get("end") or {}).get("value")
            if ftype == "signal" and e:
                sig_end = max(sig_end, int(e))
            elif ftype == "glycosylation" and s:
                desc = (f.get("description") or "").lower()
                if "n-linked" in desc or not desc:
                    glyc.append(int(s))
            elif "disulfide" in ftype and s and e:
                ss.append((int(s), int(e)))
        rows.append({"accession": acc, "signal_end": sig_end,
                     "glyco_pos": ";".join(map(str, sorted(glyc))),
                     "disulfide_pairs": ";".join(f"{a}-{b}" for a, b in ss)})
    df = pd.DataFrame(rows)
    print(f"[fetch] UniProt annotations for {len(df)} accessions; "
          f"{int((df['signal_end'] > 0).sum())} have a signal peptide")
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
