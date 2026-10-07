"""Paths and knobs. Edit here, not in the analysis modules."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = DATA / "cache"
RESULTS = ROOT / "results"
FIGS = ROOT / "figures"
for _p in (DATA, CACHE, RESULTS, FIGS):
    _p.mkdir(parents=True, exist_ok=True)

# --- data selection -------------------------------------------------------
GPCRDB = "https://gpcrdb.org/services"
CLASS_A_SLUG = "001"          # GPCRdb slug for Class A (rhodopsin-like)
SPECIES = "Homo sapiens"      # set to None for all species
SOURCE_KEEP = ("SWISSPROT",)  # drop TrEMBL entries; set to None to keep all
INCLUDE_OLFACTORY = False     # ~400 ORs; they swamp everything if included
INCLUDE_TASTE = False
# Orphans are NOT filtered: they sit inside the Class A tree and come along
# automatically with the descendants query.

USE_UNIPROT_ANNOT = True      # signal peptides, real glycosites, real disulfides
TRIM_SIGNAL_PEPTIDE = True    # the mature receptor is what a cell displays

# --- segments -------------------------------------------------------------
LOOPS = ["ECL1", "ECL2", "ECL3"]
SEGMENTS = ["N-term"] + LOOPS          # everything extracellular
PLOT_SEGMENTS = SEGMENTS               # N-term plotted everywhere the loops are
ECF = "ECF"                            # N-term+ECL1+ECL2+ECL3 concatenated
ALL_SEGMENTS = SEGMENTS + [ECF]

MIN_LOOP_LEN = 3
KS = (5, 8)                   # k-mer sizes. 8 ~ linear B-cell epitope core
ID_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
CROSSREACT_PID = 0.70
LENGTH_MARKS = [10, 25, 50, 100]       # reference lines on the length plots

# --- alignment ------------------------------------------------------------
GAP_OPEN = -11.0
GAP_EXTEND = -1.0
FREE_END_GAPS = True

# --- display tractability heuristic ---------------------------------------
PEPTIDE_LEN_OK = (8, 45)      # comfortable pIII fusion peptide range
PEPTIDE_LEN_HARD = (5, 80)    # outside this, peptide display is a bad bet
AGG_WINDOW = 7                # window for hydrophobic patch scan
FOLD_WINDOW = 21              # windowed FoldIndex for long N-termini
LRR_MIN_MOTIFS = 3            # leucine-rich repeats -> real folded ectodomain

# --- compute --------------------------------------------------------------
N_FETCH_WORKERS = 8
N_ALIGN_WORKERS = None        # None -> os.cpu_count()
REQUEST_PAUSE = 0.05

# --- plotting -------------------------------------------------------------
PALETTE = {
    "N-term": "#8172B3",
    "ECL1": "#4C72B0",
    "ECL2": "#DD8452",
    "ECL3": "#55A868",
    "ECF": "#937860",
}
FORMAT_PALETTE = {
    "linear peptide": "#4C72B0",
    "cyclic / disulfide-constrained": "#DD8452",
    "folded ectodomain": "#8172B3",
    "cell-surface only": "#C44E52",
}
