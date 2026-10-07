"""Regenerate notebooks/*.ipynb from the cell definitions below.

The notebooks are the source of truth for figures and are meant to be edited by
hand. This script exists so the set can be rebuilt if one gets mangled; it
overwrites, so do not run it casually after editing a notebook.

    python tools/build_notebooks.py
"""
from __future__ import annotations

import json
from pathlib import Path

NB = Path(__file__).resolve().parents[1] / "notebooks"
NB.mkdir(exist_ok=True)

PRELUDE = '''from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))
from src import style

style.apply()
if not style.available():
    print("Arial not found, matplotlib fell back to DejaVu Sans")

RESULTS, FIGS = ROOT / "results", ROOT / "figures"
FIGS.mkdir(exist_ok=True)
PAL = style.PALETTE
SEGMENTS = ["N-term", "ECL1", "ECL2", "ECL3"]'''

STRIP = '''def strip(ax, groups, values_by_group, log=False):
    """Jittered points with a median bar. Used in several panels."""
    for i, g in enumerate(groups):
        x = np.asarray(values_by_group[g], dtype=float)
        x = x[~np.isnan(x)]
        if not len(x):
            continue
        ax.scatter(np.random.normal(i, 0.09, len(x)), x, s=6, alpha=0.45,
                   color=PAL.get(g, "0.4"), edgecolors="none")
        ax.hlines(np.median(x), i - 0.28, i + 0.28, color="k", lw=1.6)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)
    if log:
        ax.set_yscale("log")'''

PIDHELP = '''def pid_matrix(seg):
    return pd.read_csv(RESULTS / f"pid_{seg}.csv.gz", index_col=0)


def pid_values(seg):
    v = pid_matrix(seg).values.astype(float)
    iu = np.triu_indices_from(v, k=1)
    return v[iu][~np.isnan(v[iu])]'''


# --------------------------------------------------------------------------
FIGURES = {
"01_overview": {
 "title": "Figure 1 — overview of the extracellular antigenic space",
 "intro": ("Length, pairwise similarity, how many distinct bins exist, and how "
           "much of the linear epitope space is private to a single receptor. "
           "All four extracellular segments, N-terminus included.\n\n"
           "Needs `results/` populated: run `python run_all.py` first."),
 "cells": [
  ("md", "## Load"),
  ("code", '''seg = pd.read_csv(RESULTS / "segments_long.csv")
cc = pd.read_csv(RESULTS / "cluster_counts.csv")
km = pd.read_csv(RESULTS / "kmer_per_receptor.csv")
K = int(km["k"].max())
seg.head()'''),
  ("md", "## Panels"),
  ("code", '''fig, ax = plt.subplot_mosaic([["len", "ecdf"], ["clust", "kmer"]],
                             figsize=(10.5, 7.8), constrained_layout=True)

strip(ax["len"], SEGMENTS,
      {s: seg.loc[seg.segment == s, "length"].values for s in SEGMENTS}, log=True)
ax["len"].set_ylabel("length (aa)")
ax["len"].set_title("Segment length")

a = ax["ecdf"]
for s in SEGMENTS:
    v = np.sort(pid_values(s))
    a.plot(v, np.arange(1, len(v) + 1) / len(v), color=PAL[s], label=s)
a.axvline(0.70, color="0.5", ls="--", lw=1)
a.set_xlabel("pairwise identity (shorter segment)")
a.set_ylabel("ECDF of pairs")
a.set_title("Pairwise similarity")
a.legend()

a = ax["clust"]
for s in SEGMENTS:
    d = cc[cc.segment == s].sort_values("identity_threshold")
    a.plot(d["identity_threshold"], d["n_clusters"], "o-", ms=4, color=PAL[s], label=s)
a.invert_xaxis()
a.set_xlabel("identity threshold")
a.set_ylabel("n clusters")
a.set_title("Distinct bins")
a.legend()

d = km[km.k == K]
strip(ax["kmer"], SEGMENTS,
      {s: d.loc[d.segment == s, "frac_private"].values for s in SEGMENTS})
ax["kmer"].set_ylabel(f"fraction of {K}-mers unique to that receptor")
ax["kmer"].set_title(f"Private linear epitope space (k={K})")

fig.savefig(FIGS / "01_overview.png")'''),
 ],
 "helpers": [STRIP, PIDHELP],
},

"02_targets": {
 "title": "Figure 2 — picking a target",
 "intro": ("Specificity headroom per receptor, the extracellular cysteine "
           "architecture, and glycan masking potential. The first panel is the "
           "actionable one: receptors with low identity to their nearest "
           "neighbour are the clean targets."),
 "cells": [
  ("md", "## Load"),
  ("code", '''nn = pd.read_csv(RESULTS / "nearest_neighbours.csv")
dis = pd.read_csv(RESULTS / "disulfide_audit.csv")
seg = pd.read_csv(RESULTS / "segments_long.csv")
N_LABEL = 8
nn.head()'''),
  ("md", "## Panels"),
  ("code", '''fig, ax = plt.subplot_mosaic([["nn", "nn"], ["cys", "glyc"]],
                             figsize=(10.5, 7.8), constrained_layout=True)

a = ax["nn"]
for i, s in enumerate(SEGMENTS):
    d = nn[nn.segment == s].sort_values("nearest_pid")
    y = d["nearest_pid"].values
    a.scatter(np.random.normal(i, 0.09, len(y)), y, s=7, alpha=0.5,
              color=PAL[s], edgecolors="none")
    for _, r in d.head(N_LABEL).iterrows():
        a.annotate(r["entry_name"].split("_")[0], (i + 0.16, r["nearest_pid"]),
                   fontsize=6, color="0.25")
a.set_xticks(range(len(SEGMENTS)))
a.set_xticklabels(SEGMENTS)
a.set_ylabel("identity to closest other receptor")
a.set_title("Specificity headroom (low = isolated = clean target)")

a = ax["cys"]
bars = {
    "Cys3.25 present": dis["has_C3x25"].mean(),
    "Cys in ECL2": (dis["n_cys_ECL2"] >= 1).mean(),
    "Canonical pair": dis["canonical_disulfide"].mean(),
    "2+ Cys in ECL2": dis["extra_ECL2_cys"].mean(),
    "N-term + ECL3 Cys": dis["possible_Nterm_ECL3_bond"].mean(),
}
a.barh(list(bars)[::-1], [v * 100 for v in list(bars.values())[::-1]],
       color=PAL["ECL2"], alpha=0.85)
a.set_xlabel("% of receptors")
a.set_xlim(0, 100)
a.set_title("Extracellular cysteine architecture")

a = ax["glyc"]
for s in SEGMENTS:
    d = seg[seg.segment == s]
    a.scatter(d["length"], d["n_sequon"] + np.random.normal(0, 0.06, len(d)),
              s=7, alpha=0.45, color=PAL[s], label=s, edgecolors="none")
a.set_xscale("log")
a.set_xlabel("length (aa)")
a.set_ylabel("N-X-S/T sequons")
a.set_title("Glycan masking potential")
a.legend()

fig.savefig(FIGS / "02_targets.png")'''),
 ],
 "helpers": [],
},

"03_lengths": {
 "title": "Figure 3 — are the loops really all short",
 "intro": ("Direct answer to that question. The three ECLs are peptide-scale; "
           "the N-terminus is not, and spans everything from short stubs to the "
           "LRR ectodomains of LGR4/5/6 and the glycoprotein hormone receptors. "
           "Signal peptides have already been trimmed upstream, so these are "
           "mature-receptor lengths."),
 "cells": [
  ("md", "## Load"),
  ("code", '''seg = pd.read_csv(RESULTS / "segments_long.csv")
MARKS = [10, 25, 50, 100]

seg[seg.segment.isin(SEGMENTS)].groupby("segment")["length"].describe()[
    ["count", "min", "25%", "50%", "75%", "max"]].loc[SEGMENTS]'''),
  ("md", "## Panels"),
  ("code", '''fig, ax = plt.subplot_mosaic([["hist", "hist"], ["ecdf", "frac"]],
                             figsize=(10.5, 7.5), constrained_layout=True)

a = ax["hist"]
bins = np.logspace(0, np.log10(seg["length"].max() * 1.1), 45)
for s in SEGMENTS:
    a.hist(seg.loc[seg.segment == s, "length"], bins=bins, alpha=0.5,
           color=PAL[s], label=s)
for m in MARKS:
    a.axvline(m, color="0.6", lw=0.8, ls=":")
    a.annotate(str(m), (m, a.get_ylim()[1]), fontsize=7, color="0.4",
               ha="center", va="bottom")
a.set_xscale("log")
a.set_xlabel("length (aa)")
a.set_ylabel("receptors")
a.set_title("Length distribution")
a.legend()

a = ax["ecdf"]
for s in SEGMENTS:
    x = np.sort(seg.loc[seg.segment == s, "length"].dropna().values)
    a.step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=PAL[s], label=s)
for m in MARKS:
    a.axvline(m, color="0.6", lw=0.8, ls=":")
a.set_xscale("log")
a.set_xlabel("length (aa)")
a.set_ylabel("cumulative fraction")
a.set_title("What fraction is short")
a.legend()

a = ax["frac"]
w = 0.8 / len(SEGMENTS)
for i, s in enumerate(SEGMENTS):
    x = seg.loc[seg.segment == s, "length"].dropna().values
    a.bar(np.arange(len(MARKS)) + i * w, [(x <= m).mean() * 100 for m in MARKS],
          width=w, color=PAL[s], alpha=0.85, label=s)
a.set_xticks(np.arange(len(MARKS)) + 0.4 - w / 2)
a.set_xticklabels([f"<={m}" for m in MARKS])
a.set_ylabel("% of receptors")
a.set_ylim(0, 100)
a.set_title("Fraction at or below each length")
a.legend(fontsize=7)

fig.savefig(FIGS / "03_lengths.png")'''),
 ],
 "helpers": [],
},

"04_diversity": {
 "title": "Figure 4 — how diverse are the N-termini and ECLs",
 "intro": ("Two views of diversity. The PCA asks whether the segments occupy "
           "different physicochemical territory from each other; the clustered "
           "identity heatmaps show the raw similarity structure, including "
           "whether it is driven by a few tight families against an otherwise "
           "flat background.\n\n"
           "`HEAT` picks which two segments get heatmaps."),
 "cells": [
  ("md", "## Load"),
  ("code", '''seg = pd.read_csv(RESULTS / "segments_long.csv")
FEATS = ["length", "n_cys", "n_sequon", "charge_density", "gravy",
         "f_aromatic", "f_hydrophobic", "f_gly_pro"]
HEAT = ("N-term", "ECL2")

d = seg[seg.segment.isin(SEGMENTS)].dropna(subset=FEATS)
X = d[FEATS].copy()
X["length"] = np.log10(X["length"].clip(lower=1))
xy, ev = pca(X.values)
print(f"PC1 {ev[0]:.1%}, PC2 {ev[1]:.1%}")'''),
  ("md", "## Panels"),
  ("code", '''fig, ax = plt.subplot_mosaic([["pca", "pca"], list(HEAT)],
                             figsize=(10.5, 9), constrained_layout=True)

a = ax["pca"]
for s in SEGMENTS:
    m = (d["segment"] == s).values
    a.scatter(xy[m, 0], xy[m, 1], s=10, alpha=0.55, color=PAL[s], label=s,
              edgecolors="none")
a.set_xlabel(f"PC1 ({ev[0]:.0%})")
a.set_ylabel(f"PC2 ({ev[1]:.0%})")
a.set_title("Physicochemical space (log length, charge, GRAVY, aromaticity, "
            "Cys, sequons, G/P)")
a.legend()

for s in HEAT:
    a = ax[s]
    M = np.nan_to_num(pid_matrix(s).values.astype(float), nan=0.0)
    M = (M + M.T) / 2
    np.fill_diagonal(M, 1.0)
    D = 1 - M
    np.fill_diagonal(D, 0.0)
    order = leaves_list(linkage(squareform(D, checks=False), method="average"))
    vmax = min(1.0, np.percentile(M[np.triu_indices_from(M, 1)], 99.5))
    im = a.imshow(M[np.ix_(order, order)], cmap="magma", vmin=0, vmax=vmax)
    a.set_xticks([]); a.set_yticks([])
    a.set_title(f"{s} identity ({M.shape[0]} receptors, clustered)")
    fig.colorbar(im, ax=a, fraction=0.046, label="identity")

fig.savefig(FIGS / "04_diversity.png")'''),
 ],
 "helpers": [
  PIDHELP,
  '''from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform


def pca(X):
    """Standardise, then SVD. Returns 2D scores and explained variance."""
    X = np.asarray(X, dtype=float)
    X = np.nan_to_num(X, nan=np.nanmean(X))
    X = (X - X.mean(0)) / np.where(X.std(0) == 0, 1, X.std(0))
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    return U[:, :2] * S[:2], (S ** 2 / (S ** 2).sum())[:2]''',
 ],
},

"05_displayability": {
 "title": "Figure 5 — can it be displayed on phage",
 "intro": ("FoldIndex and the charge-hydropathy plane, plus the construct-format "
           "call from `src/express.py`.\n\n"
           "The counter-intuitive read: **disorder is good**. A segment that is "
           "disordered in the intact receptor is one whose isolated peptide "
           "samples roughly the native ensemble, so a binder raised on it has a "
           "real shot on the cell. A segment predicted to fold will not "
           "reproduce that fold as a peptide."),
 "cells": [
  ("md", "## Load"),
  ("code", '''dt = pd.read_csv(RESULTS / "display_tractability.csv")
seg = pd.read_csv(RESULTS / "segments_long.csv")
m = seg.merge(dt[["entry_name", "segment", "fold_index", "display_score",
                  "recommended_format"]], on=["entry_name", "segment"])

dt.groupby(["segment", "recommended_format"]).size().unstack(fill_value=0)'''),
  ("md", "## Panels"),
  ("code", '''fig, ax = plt.subplot_mosaic([["ch", "fi"], ["fmt", "score"]],
                             figsize=(10.5, 7.8), constrained_layout=True)

a = ax["ch"]
H = (m["gravy"] + 4.5) / 9.0
R = m["charge_density"].abs()
for s in SEGMENTS:
    k = (m["segment"] == s).values
    a.scatter(H[k], R[k], s=10, alpha=0.55, color=PAL[s], label=s, edgecolors="none")
hh = np.linspace(float(H.min()), float(H.max()), 50)
a.plot(hh, 2.785 * hh - 1.151, "k--", lw=1)
a.annotate("FoldIndex = 0", (hh[-1], 2.785 * hh[-1] - 1.151), fontsize=7,
           ha="right", va="bottom", color="0.3")
a.set_xlabel("mean scaled hydrophobicity  <H>")
a.set_ylabel("mean net charge  |<R>|")
a.set_title("Charge-hydropathy: above line = disordered")
a.legend(fontsize=7)

a = ax["fi"]
for s in SEGMENTS:
    x = m[m.segment == s]
    a.scatter(x["length"], x["fold_index"], s=10, alpha=0.55, color=PAL[s],
              label=s, edgecolors="none")
a.axhline(0, color="0.5", lw=1, ls="--")
a.set_xscale("log")
a.set_xlabel("length (aa)")
a.set_ylabel("FoldIndex")
a.set_title("Below 0 = disordered = peptide is faithful")
a.legend(fontsize=7)

a = ax["fmt"]
ct = (dt.groupby(["segment", "recommended_format"]).size()
        .unstack(fill_value=0).reindex(SEGMENTS).fillna(0))
bottom = np.zeros(len(ct))
for col in ct.columns:
    a.bar(ct.index, ct[col].values, bottom=bottom,
          color=style.FORMAT_PALETTE.get(col, "0.6"), label=col)
    bottom += ct[col].values
a.set_ylabel("receptors")
a.set_title("Recommended construct format")
a.legend(fontsize=6.5)

strip(ax["score"], SEGMENTS,
      {s: dt.loc[dt.segment == s, "display_score"].values for s in SEGMENTS})
ax["score"].set_ylabel("display tractability score")
ax["score"].set_title("Composite (heuristic — see src/express.py)")

fig.savefig(FIGS / "05_displayability.png")'''),
 ],
 "helpers": [STRIP],
},
}


def cell(kind, src):
    src = src.rstrip("\n").split("\n")
    src = [l + "\n" for l in src[:-1]] + [src[-1]]
    if kind == "md":
        return {"cell_type": "markdown", "metadata": {}, "source": src}
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": src}


def build(name, spec):
    cells = [cell("md", f"# {spec['title']}\n\n{spec['intro']}"),
             cell("md", "## Setup"),
             cell("code", PRELUDE)]
    for h in spec["helpers"]:
        cells.append(cell("code", h))
    for kind, src in spec["cells"]:
        cells.append(cell(kind, src))
    nb = {"cells": cells, "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python"}},
        "nbformat": 4, "nbformat_minor": 5}
    p = NB / f"{name}.ipynb"
    p.write_text(json.dumps(nb, indent=1) + "\n")
    return p


if __name__ == "__main__":
    for name, spec in FIGURES.items():
        print("wrote", build(name, spec))
