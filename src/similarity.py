"""Step 3 - how much distinct antigenic space is there, per loop.

Three complementary views, because "different" means different things to an
antibody:

  1. pairwise sequence identity  -> overall relatedness of a loop between two
     receptors. Normalised by the SHORTER loop, because a 6-mer that is fully
     contained in a 30-mer is a real cross-reactivity liability even though a
     global alignment would score it as mostly gap.
  2. k-mer sharing (k=8)         -> linear-epitope-level overlap. A distinct
     8-mer is roughly the footprint of a linear paratope contact patch.
  3. clustering at identity cut-offs -> the headline number: how many mutually
     dissimilar loop bins exist across the whole class.
"""
from __future__ import annotations

import itertools
import os
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from . import config as C

_SEQS = None
_ALN = None


def make_aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices
    a = PairwiseAligner()
    a.mode = "global"
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = C.GAP_OPEN
    a.extend_gap_score = C.GAP_EXTEND
    if C.FREE_END_GAPS:
        try:
            a.end_gap_score = 0.0
        except AttributeError:  # older Biopython
            a.target_end_gap_score = 0.0
            a.query_end_gap_score = 0.0
    return a


def identity_count(a: str, b: str, aligner) -> int:
    if not a or not b:
        return 0
    aln = aligner.align(a, b)[0]
    ta, qa = aln.aligned
    n = 0
    for (s1, e1), (s2, e2) in zip(ta, qa):
        n += sum(1 for x, y in zip(a[s1:e1], b[s2:e2]) if x == y)
    return n


def _init(seqs):
    global _SEQS, _ALN
    _SEQS = seqs
    _ALN = make_aligner()


def _row(i):
    n = len(_SEQS)
    out = np.zeros(n, dtype=np.int32)
    a = _SEQS[i]
    for j in range(i + 1, n):
        out[j] = identity_count(a, _SEQS[j], _ALN)
    return i, out


def pid_matrix(names, seqs, workers=None):
    """Percent identity normalised by the shorter of the two loops."""
    n = len(seqs)
    ident = np.zeros((n, n), dtype=np.int32)
    workers = workers or C.N_ALIGN_WORKERS or os.cpu_count() or 1
    if workers > 1 and n > 40:
        with Pool(workers, initializer=_init, initargs=(seqs,)) as p:
            for i, row in p.imap_unordered(_row, range(n), chunksize=4):
                ident[i] = row
    else:
        _init(seqs)
        for i in range(n):
            _, ident[i] = _row(i)
    ident = np.maximum(ident, ident.T)
    L = np.array([len(s) for s in seqs], dtype=float)
    denom = np.minimum.outer(L, L)
    with np.errstate(divide="ignore", invalid="ignore"):
        pid = np.where(denom > 0, ident / denom, np.nan)
    np.fill_diagonal(pid, 1.0)
    return pd.DataFrame(pid, index=names, columns=names)


# --------------------------------------------------------------------------
def kmer_stats(names, seqs, k):
    sets = {nm: {s[i:i + k] for i in range(len(s) - k + 1)}
            for nm, s in zip(names, seqs)}
    counts = Counter()
    for v in sets.values():
        counts.update(v)
    total = len(counts)
    private = {km for km, c in counts.items() if c == 1}
    per = []
    for nm in names:
        ks = sets[nm]
        per.append({"entry_name": nm, "k": k, "n_kmer": len(ks),
                    "n_private": len(ks & private),
                    "frac_private": (len(ks & private) / len(ks)) if ks else np.nan})
    summary = {"k": k, "n_distinct_kmers": total,
               "n_private_kmers": len(private),
               "frac_private_kmers": (len(private) / total) if total else np.nan}
    return pd.DataFrame(per), summary, sets


def cluster_counts(pid: pd.DataFrame, thresholds=None):
    thresholds = thresholds or C.ID_THRESHOLDS
    D = 1.0 - pid.values.astype(float)
    D = np.nan_to_num(D, nan=1.0)
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method="average")
    rows, labels = [], {}
    for t in thresholds:
        lab = fcluster(Z, t=1.0 - t, criterion="distance")
        labels[t] = lab
        rows.append({"identity_threshold": t, "n_clusters": len(set(lab)),
                     "n_singletons": sum(v == 1 for v in Counter(lab).values())})
    return pd.DataFrame(rows), labels, Z


def nearest_neighbours(pid: pd.DataFrame):
    v = pid.values.copy().astype(float)
    np.fill_diagonal(v, -np.inf)
    idx = np.nanargmax(np.nan_to_num(v, nan=-np.inf), axis=1)
    return pd.DataFrame({
        "entry_name": pid.index,
        "nearest": pid.columns[idx],
        "nearest_pid": [v[i, j] for i, j in enumerate(idx)],
    })


def top_pairs(pid: pd.DataFrame, thresh=None, cap=500):
    thresh = C.CROSSREACT_PID if thresh is None else thresh
    v = pid.values
    iu = np.triu_indices_from(v, k=1)
    m = np.nan_to_num(v[iu], nan=0.0) >= thresh
    rows = [{"a": pid.index[i], "b": pid.columns[j], "pid": v[i, j]}
            for i, j in zip(iu[0][m], iu[1][m])]
    if not rows:
        return pd.DataFrame(columns=["a", "b", "pid"])
    return (pd.DataFrame(rows).sort_values("pid", ascending=False)
            .head(cap).reset_index(drop=True))


# --------------------------------------------------------------------------
def run(seg: pd.DataFrame = None, segments=None):
    if seg is None:
        seg = pd.read_csv(C.RESULTS / "segments_long.csv")
    segments = segments or (C.SEGMENTS + [C.ECF])

    clus_rows, nn_all, km_rows, km_summ = [], [], [], []
    for s in segments:
        sub = seg[(seg["segment"] == s) & (seg["length"] >= C.MIN_LOOP_LEN)]
        sub = sub.dropna(subset=["seq"])
        names = sub["entry_name"].tolist()
        seqs = sub["seq"].tolist()
        if len(names) < 3:
            continue
        print(f"[sim] {s}: {len(names)} receptors, "
              f"{len(names)*(len(names)-1)//2:,} pairs")
        pid = pid_matrix(names, seqs)
        pid.to_csv(C.RESULTS / f"pid_{s}.csv.gz", compression="gzip")

        cc, _, _ = cluster_counts(pid)
        cc["segment"] = s
        clus_rows.append(cc)

        nn = nearest_neighbours(pid)
        nn["segment"] = s
        nn_all.append(nn)

        top_pairs(pid).assign(segment=s).to_csv(
            C.RESULTS / f"crossreactive_pairs_{s}.csv", index=False)

        for k in C.KS:
            per, summ, _ = kmer_stats(names, seqs, k)
            per["segment"] = s
            summ["segment"] = s
            km_rows.append(per)
            km_summ.append(summ)

    pd.concat(clus_rows).to_csv(C.RESULTS / "cluster_counts.csv", index=False)
    pd.concat(nn_all).to_csv(C.RESULTS / "nearest_neighbours.csv", index=False)
    pd.concat(km_rows).to_csv(C.RESULTS / "kmer_per_receptor.csv", index=False)
    pd.DataFrame(km_summ).to_csv(C.RESULTS / "kmer_summary.csv", index=False)
    print("[sim] done")


if __name__ == "__main__":
    run()
