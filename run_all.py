"""Run the whole thing:  python run_all.py

Steps are independent and cache to disk, so you can re-run any single one:
    python -m src.fetch
    python -m src.loops
    python -m src.similarity
    python -m src.report

Figures are notebooks, one per figure, in notebooks/ - edit those directly.
"""
import argparse
import subprocess
import sys
import time

from src import config as C
from src import express, fetch, loops, report, similarity


def run_notebooks(skip=False):
    """Execute each figure notebook headlessly. The notebooks are the source of
    truth for figures; this just runs them so run_all.py produces the PNGs.
    Outputs are not written back into the .ipynb, so git stays clean."""
    if skip:
        return
    nbs = sorted((C.ROOT / "notebooks").glob("*.ipynb"))
    if not nbs:
        print("[figs] no notebooks found")
        return
    for nb in nbs:
        r = subprocess.run(
            [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook",
             "--execute", "--stdout", str(nb)],
            capture_output=True, text=True, cwd=C.ROOT)
        print(f"[figs] {nb.name}: {'ok' if r.returncode == 0 else 'FAILED'}")
        if r.returncode:
            print(r.stderr[-1500:])
    print(f"[figs] -> {C.FIGS}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-fetch", action="store_true",
                    help="reuse data/residues.csv.gz")
    ap.add_argument("--limit", type=int, default=None,
                    help="only N receptors, for a smoke test")
    ap.add_argument("--no-figures", action="store_true",
                    help="skip notebook execution")
    ap.add_argument("--no-uniprot", action="store_true",
                    help="skip UniProt annotations entirely")
    a = ap.parse_args()

    t0 = time.time()
    import pandas as pd
    annot = None
    if a.skip_fetch:
        prot = pd.read_csv(C.DATA / "proteins.csv")
        res = pd.read_csv(C.DATA / "residues.csv.gz")
        annot_path = C.DATA / "uniprot_annotations.csv"
        if annot_path.exists():
            annot = pd.read_csv(annot_path)
    else:
        prot = fetch.class_a_proteins()
        if a.limit:
            prot = prot.head(a.limit)
        print(f"[fetch] {len(prot)} Class A receptors "
              f"({int(prot['is_orphan'].sum())} flagged orphan)")
        res = fetch.residue_table(prot["entry_name"].tolist())
        if C.USE_UNIPROT_ANNOT and not a.no_uniprot:
            try:
                annot = fetch.uniprot_annotations(
                    prot["accession"].dropna().tolist())
            except Exception as e:
                print(f"[fetch] UniProt step failed ({e.__class__.__name__}); "
                      "continuing without signal-peptide trimming")

    seg, dis, wide, glyc = loops.run(res, prot, annot)
    express.run(seg, glyc)
    similarity.run(seg)
    run_notebooks(skip=a.no_figures)
    report.run()
    print(f"\ndone in {time.time()-t0:.0f}s -> {C.RESULTS}, {C.FIGS}")


if __name__ == "__main__":
    main()
