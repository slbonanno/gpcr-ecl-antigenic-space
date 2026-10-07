"""Step 5 - collapse everything into results/REPORT.md with the actual numbers."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C


def _pid_vals(seg):
    p = pd.read_csv(C.RESULTS / f"pid_{seg}.csv.gz", index_col=0)
    v = p.values.astype(float)
    iu = np.triu_indices_from(v, k=1)
    v = v[iu]
    return v[~np.isnan(v)]


def run():
    seg = pd.read_csv(C.RESULTS / "segments_long.csv")
    dis = pd.read_csv(C.RESULTS / "disulfide_audit.csv")
    cc = pd.read_csv(C.RESULTS / "cluster_counts.csv")
    kms = pd.read_csv(C.RESULTS / "kmer_summary.csv")
    nn = pd.read_csv(C.RESULTS / "nearest_neighbours.csv")
    prot = pd.read_csv(C.DATA / "proteins.csv")

    L = []
    L.append("# Class A GPCR extracellular loops - antigenic space\n")
    L.append(f"Receptors analysed: **{seg['entry_name'].nunique()}** "
             f"(orphan-family members: {int(prot['is_orphan'].sum())}; "
             f"olfactory included: {C.INCLUDE_OLFACTORY})\n")

    L.append("\n## 1. Loop geometry\n")
    L.append("| segment | n | median len | IQR | min | max | median sequons |")
    L.append("|---|---|---|---|---|---|---|")
    for s in C.SEGMENTS + [C.ECF]:
        d = seg[seg["segment"] == s]
        if not len(d):
            continue
        q1, q3 = d["length"].quantile([0.25, 0.75])
        L.append(f"| {s} | {len(d)} | {d['length'].median():.0f} | "
                 f"{q1:.0f}-{q3:.0f} | {d['length'].min():.0f} | "
                 f"{d['length'].max():.0f} | {d['n_sequon'].median():.0f} |")

    L.append("\n## 2. Pairwise similarity\n")
    L.append("Identity normalised by the shorter loop of each pair.\n")
    L.append("| segment | median pid | 90th pct | % pairs >= "
             f"{C.CROSSREACT_PID:.0%} |")
    L.append("|---|---|---|---|")
    for s in C.SEGMENTS + [C.ECF]:
        try:
            v = _pid_vals(s)
        except FileNotFoundError:
            continue
        L.append(f"| {s} | {np.median(v):.3f} | {np.percentile(v, 90):.3f} | "
                 f"{100 * (v >= C.CROSSREACT_PID).mean():.2f}% |")

    L.append("\n## 3. How many distinct bins\n")
    piv = cc.pivot(index="identity_threshold", columns="segment",
                   values="n_clusters").sort_index(ascending=False)
    L.append("| identity cut-off | " + " | ".join(map(str, piv.columns)) + " |")
    L.append("|---" * (len(piv.columns) + 1) + "|")
    for t, row in piv.iterrows():
        L.append(f"| {t:.0%} | " + " | ".join(f"{int(x)}" for x in row) + " |")

    L.append("\n## 4. Linear epitope space (k-mers)\n")
    L.append("| segment | k | distinct k-mers | unique to one receptor |")
    L.append("|---|---|---|---|")
    for _, r in kms.sort_values(["segment", "k"]).iterrows():
        L.append(f"| {r['segment']} | {int(r['k'])} | {int(r['n_distinct_kmers'])} | "
                 f"{r['frac_private_kmers']:.1%} |")

    L.append("\n## 5. Cysteines / disulfides\n")
    L.append(f"- Cys at 3.25 (top of TM3): **{100*dis['has_C3x25'].mean():.1f}%**")
    L.append(f"- >=1 Cys in ECL2: **{100*(dis['n_cys_ECL2']>=1).mean():.1f}%**")
    L.append(f"- Both, i.e. canonical TM3-ECL2 disulfide: "
             f"**{100*dis['canonical_disulfide'].mean():.1f}%**")
    L.append(f"- >=2 Cys in ECL2 (second internal bond): "
             f"{100*dis['extra_ECL2_cys'].mean():.1f}%")
    L.append(f"- Cys in both N-term and ECL3: "
             f"{100*dis['possible_Nterm_ECL3_bond'].mean():.1f}%\n")
    miss = dis[~dis["canonical_disulfide"]]["entry_name"].tolist()
    L.append(f"Receptors lacking the canonical pair ({len(miss)}): "
             + ", ".join(sorted(miss)) + "\n")

    L.append("\n## 6. Most isolated receptors per loop\n")
    L.append("Lowest identity to their nearest neighbour = cleanest specificity "
             "headroom for a loop-directed binder.\n")
    for s in C.SEGMENTS:
        d = nn[nn["segment"] == s].sort_values("nearest_pid").head(12)
        L.append(f"\n**{s}**\n")
        L.append("| receptor | nearest | pid |")
        L.append("|---|---|---|")
        for _, r in d.iterrows():
            L.append(f"| {r['entry_name']} | {r['nearest']} | "
                     f"{r['nearest_pid']:.2f} |")

    # ---- displayability -------------------------------------------------
    dt_path = C.RESULTS / "display_tractability.csv"
    if dt_path.exists():
        dt = pd.read_csv(dt_path)
        L.append("\n## 7. Peptide display tractability\n")
        L.append("Heuristic triage, weights in `src/express.py`. FoldIndex<0 "
                 "means the segment is disordered in the intact receptor, which "
                 "is what you want: an isolated peptide then samples roughly the "
                 "native ensemble.\n")
        L.append("| segment | median FoldIndex | % disordered | median score |")
        L.append("|---|---|---|---|")
        for s in C.SEGMENTS:
            d = dt[dt["segment"] == s]
            if not len(d):
                continue
            L.append(f"| {s} | {d['fold_index'].median():.3f} | "
                     f"{100*(d['fold_index']<0).mean():.0f}% | "
                     f"{d['display_score'].median():.2f} |")
        L.append("\n**Recommended construct format**\n")
        ct = (dt.groupby(["segment", "recommended_format"]).size()
                .unstack(fill_value=0))
        L.append("| segment | " + " | ".join(ct.columns) + " |")
        L.append("|---" * (len(ct.columns) + 1) + "|")
        for s, row in ct.iterrows():
            L.append(f"| {s} | " + " | ".join(str(int(x)) for x in row) + " |")
        for s in C.SEGMENTS:
            d = dt[dt["segment"] == s].nlargest(12, "display_score")
            if not len(d):
                continue
            L.append(f"\n**Best {s} peptide-display candidates**\n")
            L.append("| receptor | len | FoldIndex | score | format | flags |")
            L.append("|---|---|---|---|---|---|")
            for _, r in d.iterrows():
                L.append(f"| {r['entry_name']} | {int(r['length'])} | "
                         f"{r['fold_index']:.2f} | {r['display_score']:.2f} | "
                         f"{r['recommended_format']} | {r['flags'] if isinstance(r['flags'],str) else ''} |")

    (C.RESULTS / "REPORT.md").write_text("\n".join(L))
    print(f"[report] {C.RESULTS/'REPORT.md'}")


if __name__ == "__main__":
    run()
