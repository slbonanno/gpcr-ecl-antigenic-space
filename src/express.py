"""Can this extracellular segment be shown to a phage library as a peptide,
and will binders selected on it still bind the receptor on a cell?

This is a transparent triage heuristic, not a predictor. Every sub-score is a
published or mechanistically obvious quantity and all of them are written out
individually, so you can disagree with the weighting and re-rank yourself.

The reasoning, in order of how much it actually matters:

FOLD STATE.  FoldIndex (Prilusky 2005) = 2.785*<H> - |<R>| - 1.151, on the
  Kyte-Doolittle scale rescaled to 0-1. Negative = intrinsically disordered.
  For peptide display the counter-intuitive result is that DISORDER IS GOOD:
  if the segment is disordered in the intact receptor, a synthetic peptide of
  it samples roughly the same conformational ensemble as the native thing, so
  a binder raised on the peptide has a real chance on the cell. If the segment
  is predicted to fold, the isolated peptide will not reproduce that fold and
  you will select against a conformation that does not exist in vivo.

LENGTH.  Below ~8 aa there is not enough surface for a paratope to grip and
  what you get is dominated by flanking vector sequence. Above ~45-80 aa you
  are no longer displaying a peptide, you are displaying a domain, and it needs
  to fold correctly in the periplasm to be worth anything.

CYSTEINE.  An odd number of free cysteines in an oxidising periplasm scrambles.
  An even number is an opportunity: cyclise it and you recover the native
  constraint. Zero is the easiest to make and the least native-like.

GLYCANS.  Phage display in E. coli is aglycosylated. A segment carrying several
  N-glycans is, on a real cell, substantially covered. Binders selected against
  the bare backbone can be selecting epitopes that are physically inaccessible
  on the receptor. This is the most common reason a beautiful peptide-derived
  clone does nothing in a whole-cell binding assay.

AGGREGATION.  Contiguous hydrophobic, low-net-charge windows drive fusion
  toxicity, poor phage yield, and enrichment of amyloid-binders.

LEUCINE-RICH REPEATS.  LGR4/5/6 and the glycoprotein hormone receptors carry
  genuine LRR ectodomains. Those are proteins, not peptides. Detected and
  routed to a different construct recommendation rather than penalised.
"""
from __future__ import annotations

import re
import numpy as np
import pandas as pd

from . import config as C
from .loops import KD, SEQUON

LRR = re.compile(r"L..L.L..[NCTSD]")
_KDMIN, _KDRANGE = -4.5, 9.0


def _scaled_h(s):
    return np.array([(KD[a] + 4.5) / _KDRANGE for a in s if a in KD])


def _charge_vec(s):
    return np.array([1.0 if a in "KR" else (-1.0 if a in "DE" else 0.0)
                     for a in s if a in KD])


def fold_index(s: str):
    h, q = _scaled_h(s), _charge_vec(s)
    if len(h) == 0:
        return np.nan
    return 2.785 * h.mean() - abs(q.mean()) - 1.151


def windowed_fold_index(s: str, w=None):
    """Most-folded and most-disordered stretch. Useful for long N-termini
    that are part globular domain, part flexible linker."""
    w = w or C.FOLD_WINDOW
    if len(s) < w:
        fi = fold_index(s)
        return fi, fi
    vals = [fold_index(s[i:i + w]) for i in range(len(s) - w + 1)]
    return float(np.nanmax(vals)), float(np.nanmin(vals))


def agg_hotspots(s: str, w=None):
    """Max mean hydrophobicity over a low-charge window, and how many such
    windows there are."""
    w = w or C.AGG_WINDOW
    if len(s) < w:
        return np.nan, 0
    h, q = _scaled_h(s), _charge_vec(s)
    if len(h) < w:
        return np.nan, 0
    best, n = -np.inf, 0
    for i in range(len(h) - w + 1):
        hw, qw = h[i:i + w].mean(), abs(q[i:i + w].sum())
        if qw <= 1:
            best = max(best, hw)
            if hw >= 0.55:
                n += 1
    return (best if np.isfinite(best) else np.nan), n


def _ramp(x, lo, hi):
    """0 below lo, 1 above hi, linear between."""
    if hi == lo:
        return float(x >= hi)
    return float(np.clip((x - lo) / (hi - lo), 0, 1))


def score_segment(seq: str, n_glyc: int = None) -> dict:
    seq = seq or ""
    n = len(seq)
    if n == 0:
        return {}
    fi = fold_index(seq)
    fi_max, fi_min = windowed_fold_index(seq)
    agg, n_agg = agg_hotspots(seq)
    n_cys = seq.count("C")
    n_seq = len(SEQUON.findall(seq)) if n_glyc is None else n_glyc
    n_lrr = len(LRR.findall(seq))

    lo_ok, hi_ok = C.PEPTIDE_LEN_OK
    lo_hard, hi_hard = C.PEPTIDE_LEN_HARD
    s_len = min(_ramp(n, lo_hard, lo_ok), 1 - _ramp(n, hi_ok, hi_hard))

    # disorder is GOOD for faithful peptide display (see module docstring)
    s_disorder = float(np.clip(0.5 - fi * 2.0, 0, 1))
    s_cys = {0: 0.85}.get(n_cys, 1.0 if n_cys % 2 == 0 else 0.25)
    s_glyc = float(np.clip(1 - (n_seq / max(n, 1)) * 25, 0, 1))
    s_agg = float(np.clip(1 - max(n_agg - 1, 0) / 4, 0, 1))

    score = float(np.mean([s_len, s_disorder, s_cys, s_glyc, s_agg]))

    if n_lrr >= C.LRR_MIN_MOTIFS or (n > hi_hard and fi > 0):
        fmt = "folded ectodomain"
    elif n < lo_hard or s_glyc < 0.35:
        fmt = "cell-surface only"
    elif n_cys >= 2 and n_cys % 2 == 0:
        fmt = "cyclic / disulfide-constrained"
    elif n_cys % 2 == 1 or n > hi_ok:
        fmt = "cell-surface only" if n_cys % 2 == 1 and n_cys > 1 else "linear peptide"
    else:
        fmt = "linear peptide"

    flags = []
    if n_cys % 2 == 1:
        flags.append("odd Cys (scrambling risk)")
    if n_seq >= 2:
        flags.append(f"{n_seq} glycosites (masked on cell)")
    if n_agg >= 2:
        flags.append("hydrophobic patches")
    if fi > 0.05:
        flags.append("predicted folded; peptide will not reproduce it")
    if n < lo_ok:
        flags.append("too short to grip")

    return {
        "length": n, "fold_index": fi, "fi_window_max": fi_max,
        "fi_window_min": fi_min, "agg_max": agg, "n_agg_windows": n_agg,
        "n_cys": n_cys, "n_glyco": n_seq, "n_lrr_motifs": n_lrr,
        "s_length": s_len, "s_disorder": s_disorder, "s_cys": s_cys,
        "s_glycan": s_glyc, "s_aggregation": s_agg,
        "display_score": score, "recommended_format": fmt,
        "flags": "; ".join(flags),
    }


def run(seg: pd.DataFrame = None, glyc: pd.DataFrame = None):
    if seg is None:
        seg = pd.read_csv(C.RESULTS / "segments_long.csv")
    gmap = {}
    if glyc is not None and len(glyc):
        gmap = {(r.entry_name, r.segment): int(r.n_glyco_annotated)
                for r in glyc.itertuples()}

    rows = []
    for r in seg.itertuples():
        if r.segment not in C.SEGMENTS:
            continue
        d = score_segment(r.seq, gmap.get((r.entry_name, r.segment)))
        if not d:
            continue
        d["entry_name"] = r.entry_name
        d["segment"] = r.segment
        rows.append(d)
    df = pd.DataFrame(rows)
    cols = ["entry_name", "segment"] + [c for c in df.columns
                                        if c not in ("entry_name", "segment")]
    df = df[cols].sort_values(["segment", "display_score"], ascending=[True, False])
    df.to_csv(C.RESULTS / "display_tractability.csv", index=False)
    print("[express] recommended formats:")
    print(df.groupby(["segment", "recommended_format"]).size()
            .unstack(fill_value=0).to_string())
    return df


if __name__ == "__main__":
    run()
