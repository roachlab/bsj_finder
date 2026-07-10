
import sys
import logging
from pathlib import Path
#import pandas as pd

# Import detect_bsj_paf from paf_parser.py in the same scripts/ folder ─────
# __file__ = workflow/scripts/detect_bsj_tsv.py
# paf_parser.py lives in the same directory
sys.path.insert(0, str(Path(__file__).parent))##todo confirm from mike
from paf_parser import detect_bsj_paf

# Redirect logging to Snakemake log file ────────────────────────────────────
log_file = str(snakemake.log[0])
logging.basicConfig(
    filename = log_file,
    filemode = "w",
    level    = logging.INFO,
    format   = "%(levelname)s: %(message)s"
)

# Get inputs and params from Snakemake ─────────────────────────────────────
paf_file = str(snakemake.input.paf)
tsv_out  = str(snakemake.output.tsv)

# Tunable thresholds — set in rule params, fall back to defaults
min_align_len    = snakemake.params.get("min_align_len",    50)
max_bsj_gap      = snakemake.params.get("max_bsj_gap",      -3)
max_bsj_ovlp     = snakemake.params.get("max_bsj_ovlp",     15)
min_pident       = snakemake.params.get("min_pident",        0.95)
min_mapq         = snakemake.params.get("min_mapq",          20)
valid_coord_dist = snakemake.params.get("valid_coord_dist",  10)
max_gap_size     = snakemake.params.get("max_gap_size",      100000)

logging.info(f"Input PAF  : {paf_file}")
logging.info(f"Output TSV : {tsv_out}")
logging.info(f"Parameters : min_pident={min_pident}, min_mapq={min_mapq}, "
             f"max_bsj_ovlp={max_bsj_ovlp}")

# Ensure output directory exists
Path(tsv_out).parent.mkdir(parents=True, exist_ok=True)

# Run detect_bsj_paf() and write output TSV ────────────────────────────────
n_bsj = 0
with open(tsv_out, "w") as out_tsv:
    out_tsv.write(
        "q_id\tbsj_coord\tchromosome\t"
        "s_start\ts_end\ts_strand\toverlap_trimmed_nt\n"
    )
    for line in detect_bsj_paf(
        file_paf_gz      = paf_file,
        min_align_len    = min_align_len,
        max_bsj_gap      = max_bsj_gap,
        max_bsj_ovlp     = max_bsj_ovlp,
        min_pident       = min_pident,
        min_m_quality    = min_mapq,
        valid_coord_dist = valid_coord_dist,
        max_gap_size     = max_gap_size,
    ):
        out_tsv.write(line)
        n_bsj += 1

logging.info(f"BSJ entries written: {n_bsj:,} → {tsv_out}")

if n_bsj == 0:
    logging.warning(
        "No BSJ detected — check PAF file and filter parameters"
    )