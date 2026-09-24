
import sys
import gzip
import logging
from pathlib import Path

# Import the junction-sequence helpers from bsj_seq_utils.py (same scripts/ dir) ──
sys.path.insert(0, str(Path(__file__).parent))

# from workflow.scripts import extract_bsj_sequences
from bsj_seq_utils import load_bsj_coordinates, extract_bsj_sequences

# Redirect logging to Snakemake log file ────────────────────────────────────
log_file = str(snakemake.log[0])
logging.basicConfig(
    filename = log_file,
    filemode = "w",
    level    = logging.INFO,
    format   = "%(levelname)s: %(message)s",
)

# Get inputs, outputs and params from Snakemake ────────────────────────────
tsv_file     = str(snakemake.input.tsv)
r1_file      = str(snakemake.input.r1)
fasta_out    = str(snakemake.output.fasta)

r2_file      = snakemake.params.get("r2_in",      "")
orphan_file  = snakemake.params.get("orphan_in",  "")
flank        = snakemake.params.get("flank",        50)
#read_suffix  = snakemake.params.get("read_suffix",  "")

logging.info(f"Input TSV    : {tsv_file}")
logging.info(f"R1 renamed FASTQ  : {r1_file}")
logging.info(f"R2 renamed FASTQ  : {r2_file or 'not provided (SE sample)'}")
logging.info(f"Orphan renamed FASTQ  : {orphan_file or 'not provided'}")
logging.info(f"R1 FASTA output  : {fasta_out}")
logging.info(f"R2 FASTA output  : {fasta_out     or 'not applicable'}")
logging.info(f"Orphan FASTA out : {fasta_out or 'not applicable'}")
logging.info(f"Flank        : {flank} bp")
# logging.info(f"Read suffix  : {read_suffix}")

# Ensure output directory exists ───────────────────────────────────────────────
Path(fasta_out).parent.mkdir(parents=True, exist_ok=True)

bsj_coords = load_bsj_coordinates(tsv_file)
fastq_sources = [r1_file, r2_file, orphan_file]

file_pairs = [(r1_file, fasta_out)]   # R1 always included
for fq_in in [r2_file, orphan_file ]:
    if fq_in and Path(fq_in).exists():
        file_pairs.append((fq_in, fasta_out)),
        logging.info(f"Will process : {fq_in}")
    elif fq_in:
        logging.info(f"Skipping     : {fq_in} (not found — SE sample or no orphans)")
# Run extract_bsj_sequences() ───────────────────────────────────────────────
# extract_bsj_sequences(
#     bsj_tsv = tsv_file,
#     fastq_file = r1_file,
#     r2_file = r2_file,
#     orphan_file = orphan_file,
#     output_fasta   = fasta_out,
#     flank          = flank,
# )
# read_id_suffix = read_suffix,
logging.info("extract_bsj_fasta.py complete ✓")

total_extracted = 0
with gzip.open(fasta_out, "wt") as out_fasta:
    for fastq_file in fastq_sources:
        for line in extract_bsj_sequences(
            bsj_coords = bsj_coords,
            fastq_file = fastq_file,
            flank = flank,
        ):
                out_fasta.write(line)
                total_extracted += 1

logging.info(f"Total BSJ sequences extracted: {total_extracted:,} → {fasta_out}")
logging.info("extract_bsj_seq.py complete")
