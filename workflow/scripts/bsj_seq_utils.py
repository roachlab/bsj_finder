#!/usr/bin/env python3
"""
Junction-sequence helpers for extract_bsj_seq.py.

HISTORY: this file used to be the PAF-based detector -- detect_bsj_paf(), a
standalone argparse CLI and a main() driving both. That path was replaced by
sa_detect.py (SAM/SA-tag based), nothing produced the .paf it consumed, and the
code was removed, and the file was renamed from paf_parser.py (2026-09-24) since it no
longer parses anything PAF-shaped.

Only the two functions below are live:
  load_bsj_coordinates(tsv)  -- read_id -> bsj_coord, from a bsj_cor table
  extract_bsj_sequences(...) -- yield the +/-flank window around each read's junction

NOTE: no logging.basicConfig here, deliberately. This module used to configure the
root logger to stderr at import time, which silently won the race against the
importing script's own basicConfig(filename=snakemake.log[0]) -- so bsj_cor_fasta
wrote 0-byte log files. Logging setup belongs to the entry-point script.
"""
import gzip
import logging
from pathlib import Path

from Bio import SeqIO


def load_bsj_coordinates(bsj_tsv):
    opener = gzip.open(bsj_tsv, "rt") if bsj_tsv.endswith(".gz") else open(bsj_tsv)
    bsj_coords = {}
    with opener as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            bsj_coords[parts[0]] = int(parts[1])
    logging.info(f"Loaded {len(bsj_coords):,} BSJ coordinates from {bsj_tsv}")
    return bsj_coords


# Sequence extraction ───────────────────────────────────────────────────────
def extract_bsj_sequences(
    fastq_file=None,
    bsj_coords=None,
    flank=50,
):
    if not Path(fastq_file).exists():
        logging.info(f"FASTQ not present, skipping: {fastq_file}")
        return


    extracted = 0
    with gzip.open(fastq_file, "rt") as in_fq:
        for record in SeqIO.parse(in_fq, "fastq"):
            read_id = record.id
            if read_id not in bsj_coords:
                continue
            bsj_coord = bsj_coords[read_id]
            start = max(0, bsj_coord - flank)
            end   = min(len(record.seq), bsj_coord + flank)
            yield (
                f">{read_id}_BSJ_{bsj_coord}_region_{start}-{end}\n"
                f"{record.seq[start:end]}\n"
            )
            extracted += 1

    not_found = len(bsj_coords) - extracted
    logging.info(f"Extracted {extracted:,} sequences")
    if not_found > 0:
        logging.warning(f"{not_found} read IDs not found in FASTQ")
