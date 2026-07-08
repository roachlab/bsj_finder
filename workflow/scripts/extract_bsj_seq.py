
import sys
import logging
from pathlib import Path

# Import extract_bsj_sequences from paf_parser.py in same scripts/ folder ──
sys.path.insert(0, str(Path(__file__).parent))
from paf_parser import extract_bsj_sequences

# Redirect logging to Snakemake log file ────────────────────────────────────
log_file = str(snakemake.log[0])
logging.basicConfig(
    filename = log_file,
    filemode = "w",
    level    = logging.INFO,
    format   = "%(levelname)s: %(message)s"
)

# Get inputs, outputs and params from Snakemake ────────────────────────────
tsv_file     = str(snakemake.input.tsv)
fastq_file   = str(snakemake.input.fastq)
fasta_out    = str(snakemake.output.fasta)

flank        = snakemake.params.get("flank",        50)
read_suffix  = snakemake.params.get("read_suffix",  "/1")

logging.info(f"Input TSV    : {tsv_file}")
logging.info(f"Input FASTQ  : {fastq_file}")
logging.info(f"Output FASTA : {fasta_out}")
logging.info(f"Flank        : {flank} bp")
logging.info(f"Read suffix  : {read_suffix}")

# Ensure output directory exists ───────────────────────────────────────────────
Path(fasta_out).parent.mkdir(parents=True, exist_ok=True)

# Run extract_bsj_sequences() ───────────────────────────────────────────────
extract_bsj_sequences(
    bsj_tsv        = tsv_file,
    fastq_file     = fastq_file,
    output_fasta   = fasta_out,
    flank          = flank,
    read_id_suffix = read_suffix,
)

logging.info("extract_bsj_fasta.py complete ✓")