"""Offline end-to-end check on synthetic receptors (no network).

    python -m tests.smoke_test

Writes into data/ and results/ exactly like the real run, so delete those
folders afterwards (or just re-run run_all.py, which overwrites).
"""
import random

import numpy as np
import pandas as pd

from src import config as C
from src import express, loops, report, similarity

AA = "ACDEFGHIKLMNPQRSTVWY"
SEG_ORDER = ["N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4",
             "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "C-term"]
LENS = {"N-term": (8, 40), "TM1": (25, 30), "ICL1": (5, 10), "TM2": (25, 30),
        "ECL1": (5, 20), "TM3": (28, 34), "ICL2": (8, 16), "TM4": (24, 30),
        "ECL2": (10, 45), "TM5": (26, 34), "ICL3": (10, 60), "TM6": (28, 34),
        "ECL3": (5, 25), "TM7": (24, 30), "H8": (8, 14), "C-term": (10, 60)}


def synth(n=60, seed=0):
    rng = random.Random(seed)
    rows, prot = [], []
    families = ["Peptide", "Aminergic", "Lipid", "Orphan receptors", "Nucleotide"]
    for i in range(n):
        en = f"rcp{i:03d}_human"
        fam = families[i % len(families)]
        pos = 0
        for s in SEG_ORDER:
            lo, hi = LENS[s]
            L = rng.randint(lo, hi)
            seq = "".join(rng.choice(AA) for _ in range(L))
            if s == "ECL2" and rng.random() < 0.92:
                seq = "C" + seq[1:]
            for j, aa in enumerate(seq):
                pos += 1
                gn = ""
                if s == "TM3" and j == 0:
                    gn = "3.25"
                    aa = "C" if rng.random() < 0.9 else rng.choice("SFY")
                rows.append((en, pos, aa, s, gn))
        prot.append({"entry_name": en, "accession": f"P{i:05d}",
                     "name": en, "species": "Homo sapiens", "source": "SWISSPROT",
                     "family_slug": "001_001", "family_path": f"Class A | {fam}",
                     "sequence": "", "is_olfactory": False, "is_taste": False,
                     "is_orphan": fam == "Orphan receptors"})
    res = pd.DataFrame(rows, columns=["entry_name", "pos", "aa", "segment", "gn"])
    return res, pd.DataFrame(prot)


def main():
    res, prot = synth()
    prot.to_csv(C.DATA / "proteins.csv", index=False)
    res.to_csv(C.DATA / "residues.csv.gz", index=False, compression="gzip")

    seg, dis, wide, glyc = loops.run(res, prot, None)
    dt = express.run(seg, glyc)
    assert dt['display_score'].between(0, 1).all()
    assert set(dt['segment']) <= set(C.SEGMENTS)
    assert set(C.LOOPS).issubset(set(seg["segment"])), "missing loop segments"
    assert (seg.loc[seg.segment == C.ECF, "length"].values > 0).all()

    similarity.run(seg)
    pid = pd.read_csv(C.RESULTS / "pid_ECL2.csv.gz", index_col=0)
    v = pid.values
    assert np.allclose(np.diag(v), 1.0), "diagonal should be 1"
    assert np.allclose(v, v.T, equal_nan=True), "matrix not symmetric"
    assert np.nanmax(v[np.triu_indices_from(v, 1)]) <= 1.0

    report.run()
    print("\nSMOKE TEST PASSED")
    print((C.RESULTS / "REPORT.md").read_text()[:1200])


if __name__ == "__main__":
    main()
