"""Step 2 - turn the residue table into per-receptor extracellular segment
sequences plus features.

Handles the two things that quietly wreck GPCR extracellular analyses:

  1. SIGNAL PEPTIDES. UniProt (and therefore GPCRdb) numbering starts at
     residue 1 of the precursor, so the annotated N-term of many Class A
     receptors begins with a signal peptide that is cleaved and is NOT on the
     mature receptor at the cell surface. Displaying it on phage would be
     selecting binders against a sequence that does not exist in vivo. It is
     trimmed here and reported separately.
  2. The TM3(Cys3.25) <-> ECL2 disulfide, the one near-universal piece of
     extracellular architecture in Class A. Measured, not assumed.
"""
from __future__ import annotations

import re
import numpy as np
import pandas as pd

from . import config as C

SEQUON = re.compile(r"N[^P][ST]")
KD = {  # Kyte-Doolittle
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
AROM, HYDRO = set("FWY"), set("AVILMFWC")


def _charge(s):
    return sum(s.count(a) for a in "KR") + 0.1 * s.count("H") \
        - sum(s.count(a) for a in "DE")


def features(seq: str) -> dict:
    s = "".join(ch for ch in seq if ch in KD)
    n = len(s)
    if n == 0:
        return {"length": 0, "n_cys": 0, "n_sequon": 0, "net_charge": np.nan,
                "charge_density": np.nan, "gravy": np.nan, "f_aromatic": np.nan,
                "f_hydrophobic": np.nan, "f_gly_pro": np.nan}
    return {
        "length": len(seq),
        "n_cys": s.count("C"),
        "n_sequon": len(SEQUON.findall(seq)),
        "net_charge": _charge(s),
        "charge_density": _charge(s) / n,
        "gravy": sum(KD[a] for a in s) / n,
        "f_aromatic": sum(a in AROM for a in s) / n,
        "f_hydrophobic": sum(a in HYDRO for a in s) / n,
        "f_gly_pro": sum(a in "GP" for a in s) / n,
    }


def apply_signal_peptides(res: pd.DataFrame, prot: pd.DataFrame,
                          annot: pd.DataFrame) -> tuple:
    """Relabel precursor signal-peptide residues so they drop out of N-term."""
    if annot is None or not len(annot) or not C.TRIM_SIGNAL_PEPTIDE:
        return res, pd.DataFrame(columns=["entry_name", "signal_len"])
    acc = prot.set_index("entry_name")["accession"]
    sig = annot.set_index("accession")["signal_end"].to_dict()
    res = res.copy()
    res["_sig"] = res["entry_name"].map(lambda e: sig.get(acc.get(e), 0)).fillna(0)
    mask = (res["pos"] <= res["_sig"]) & (res["segment"] == "N-term")
    res.loc[mask, "segment"] = "Signal"
    tab = (res[res["segment"] == "Signal"].groupby("entry_name").size()
           .rename("signal_len").reset_index())
    print(f"[loops] signal peptide trimmed from {len(tab)} receptors "
          f"(median {tab['signal_len'].median() if len(tab) else 0:.0f} aa)")
    return res.drop(columns="_sig"), tab


def annotated_glycosites(res: pd.DataFrame, prot: pd.DataFrame,
                         annot: pd.DataFrame) -> pd.DataFrame:
    """Count experimentally/homology-annotated N-glycans per segment."""
    if annot is None or not len(annot):
        return pd.DataFrame(columns=["entry_name", "segment", "n_glyco_annotated"])
    acc = prot.set_index("entry_name")["accession"].to_dict()
    gp = {}
    for r in annot.itertuples():
        if isinstance(r.glyco_pos, str) and r.glyco_pos:
            gp[r.accession] = {int(x) for x in r.glyco_pos.split(";")}
    rows = []
    for en, sub in res.groupby("entry_name"):
        pos = gp.get(acc.get(en), set())
        if not pos:
            continue
        hit = sub[sub["pos"].isin(pos)]
        for s, d in hit.groupby("segment"):
            if s in C.SEGMENTS:
                rows.append({"entry_name": en, "segment": s,
                             "n_glyco_annotated": len(d)})
    return pd.DataFrame(rows)


def build_segments(res: pd.DataFrame) -> pd.DataFrame:
    res = res.sort_values(["entry_name", "pos"])
    keep = res[res["segment"].isin(C.SEGMENTS)]
    g = keep.groupby(["entry_name", "segment"], sort=False)
    seg = g.agg(seq=("aa", lambda x: "".join(x)),
                start=("pos", "min"),
                end=("pos", "max")).reset_index()

    order = {s: i for i, s in enumerate(C.SEGMENTS)}
    ecf = (seg.assign(_o=seg["segment"].map(order))
              .sort_values(["entry_name", "_o"])
              .groupby("entry_name")["seq"].apply("".join).reset_index())
    ecf["segment"] = C.ECF
    ecf["start"] = np.nan
    ecf["end"] = np.nan

    seg = pd.concat([seg, ecf], ignore_index=True)
    feat = pd.DataFrame([features(s) for s in seg["seq"]])
    seg = pd.concat([seg.reset_index(drop=True), feat], axis=1)
    seg.to_csv(C.RESULTS / "segments_long.csv", index=False)
    return seg


def disulfide_audit(res: pd.DataFrame) -> pd.DataFrame:
    gn325 = res[res["gn"] == "3.25"].set_index("entry_name")["aa"]
    rows = []
    for en, sub in res.groupby("entry_name"):
        seg = {k: "".join(v["aa"]) for k, v in
               sub.sort_values("pos").groupby("segment")}
        aa325 = gn325.get(en, None)
        ecl2 = seg.get("ECL2", "")
        rows.append({
            "entry_name": en,
            "res_3x25": aa325,
            "has_C3x25": aa325 == "C",
            "n_cys_ECL2": ecl2.count("C"),
            "n_cys_ECL1": seg.get("ECL1", "").count("C"),
            "n_cys_ECL3": seg.get("ECL3", "").count("C"),
            "n_cys_Nterm": seg.get("N-term", "").count("C"),
        })
    d = pd.DataFrame(rows)
    d["canonical_disulfide"] = d["has_C3x25"] & (d["n_cys_ECL2"] >= 1)
    d["extra_ECL2_cys"] = d["n_cys_ECL2"] >= 2
    d["possible_Nterm_ECL3_bond"] = (d["n_cys_Nterm"] >= 1) & (d["n_cys_ECL3"] >= 1)
    d.to_csv(C.RESULTS / "disulfide_audit.csv", index=False)
    return d


def wide_table(seg, prot, dis, sig) -> pd.DataFrame:
    piv = seg.pivot(index="entry_name", columns="segment", values="seq")
    lens = seg.pivot(index="entry_name", columns="segment", values="length")
    lens.columns = [f"len_{c}" for c in lens.columns]
    sq = seg.pivot(index="entry_name", columns="segment", values="n_sequon")
    sq.columns = [f"sequon_{c}" for c in sq.columns]
    w = (piv.join(lens).join(sq)
            .merge(dis.set_index("entry_name"), left_index=True,
                   right_index=True, how="left")
            .merge(prot.set_index("entry_name")[
                       ["accession", "name", "family_path", "is_orphan"]],
                   left_index=True, right_index=True, how="left"))
    if len(sig):
        w = w.merge(sig.set_index("entry_name"), left_index=True,
                    right_index=True, how="left")
    w = w.reset_index()
    w.to_csv(C.RESULTS / "receptors_wide.csv", index=False)
    return w


def run(res: pd.DataFrame = None, prot: pd.DataFrame = None,
        annot: pd.DataFrame = None):
    if res is None:
        res = pd.read_csv(C.DATA / "residues.csv.gz")
    if prot is None:
        prot = pd.read_csv(C.DATA / "proteins.csv")
    if annot is None and (C.DATA / "uniprot_annotations.csv").exists():
        annot = pd.read_csv(C.DATA / "uniprot_annotations.csv")

    res, sig = apply_signal_peptides(res, prot, annot)
    glyc = annotated_glycosites(res, prot, annot)
    glyc.to_csv(C.RESULTS / "glycosites_annotated.csv", index=False)

    seg = build_segments(res)
    dis = disulfide_audit(res)
    w = wide_table(seg, prot, dis, sig)
    print(f"[loops] {w.shape[0]} receptors; canonical TM3-ECL2 disulfide in "
          f"{100 * dis['canonical_disulfide'].mean():.1f}%")
    return seg, dis, w, glyc


if __name__ == "__main__":
    run()
