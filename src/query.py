"""Target-centric lookup:  python -m src.query drd2_human

For one receptor, prints per-loop sequence, the receptors whose loops look most
like it (off-target risk if you pan on that loop), and the stretches that are
unique to it in the whole Class A set (where to aim).
"""
from __future__ import annotations

import sys
from collections import Counter

import pandas as pd

from . import config as C


def unique_windows(target: str, seg: pd.DataFrame, segment: str, k: int):
    sub = seg[(seg["segment"] == segment) & (seg["length"] >= k)]
    tgt = sub.loc[sub["entry_name"] == target, "seq"]
    if tgt.empty:
        return []
    s = tgt.iloc[0]
    counts = Counter()
    for other in sub["seq"]:
        counts.update({other[i:i + k] for i in range(len(other) - k + 1)})
    out, run = [], []
    for i in range(len(s) - k + 1):
        km = s[i:i + k]
        if counts[km] == 1:
            run.append((i, km))
        else:
            if run:
                out.append((run[0][0] + 1, s[run[0][0]:run[-1][0] + k]))
                run = []
    if run:
        out.append((run[0][0] + 1, s[run[0][0]:run[-1][0] + k]))
    return out


def run(target: str, top=8, k=None):
    k = k or max(C.KS)
    seg = pd.read_csv(C.RESULTS / "segments_long.csv")
    dis = pd.read_csv(C.RESULTS / "disulfide_audit.csv").set_index("entry_name")

    if target not in set(seg["entry_name"]):
        cand = [e for e in seg["entry_name"].unique() if target.lower() in e.lower()]
        print(f"'{target}' not found. Close matches: {cand[:10]}")
        return

    print(f"\n=== {target} ===")
    if target in dis.index:
        r = dis.loc[target]
        print(f"3.25 residue: {r['res_3x25']}   ECL2 Cys: {int(r['n_cys_ECL2'])}   "
              f"canonical disulfide: {bool(r['canonical_disulfide'])}")

    for s in C.LOOPS:
        row = seg[(seg["segment"] == s) & (seg["entry_name"] == target)]
        if row.empty:
            continue
        row = row.iloc[0]
        print(f"\n--- {s}  ({int(row['length'])} aa, {int(row['n_sequon'])} sequons) ---")
        print(row["seq"])
        try:
            pid = pd.read_csv(C.RESULTS / f"pid_{s}.csv.gz", index_col=0)
            if target in pid.index:
                v = pid.loc[target].drop(target).sort_values(ascending=False)
                print("closest others:", ", ".join(
                    f"{n}({x:.2f})" for n, x in v.head(top).items()))
        except FileNotFoundError:
            pass
        uw = unique_windows(target, seg, s, k)
        if uw:
            print(f"stretches unique in Class A (>={k}aa):")
            for pos, sq in uw:
                print(f"   +{pos:>3}  {sq}")
        else:
            print(f"   no {k}-mer unique to this receptor in this loop")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
    else:
        run(sys.argv[1])
