#!/usr/bin/env python3
#Imports─────────────────────────────────────────────────────────────
import sys
import gzip
import argparse
import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO
# from workflow.scripts import extract_bsj_seq

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    stream=sys.stderr
)


# Argument parser ─────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Detect BSJ coordinates from PAF and extract sequences"
    )
    p.add_argument("--paf",   required=True)
    p.add_argument("--fastq", required=True)
    p.add_argument("--tsv",   required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--flank",        type=int,   default=50)
    p.add_argument("--read-suffix",  type=str,   default="/1")
    p.add_argument("--min-align-len", type=int, default=30)
    p.add_argument("--min-align-frac", type=float, default=0.35)
    p.add_argument("--min-align-len",type=int,   default=50)
    p.add_argument("--max-bsj-gap",  type=int,   default=-3)
    p.add_argument("--max-bsj-ovlp", type=int,   default=15)
    p.add_argument("--min-pident",   type=float, default=0.95)
    p.add_argument("--min-mapq",     type=int,   default=20)
    p.add_argument("--valid-coord-dist", type=int, default=10)
    p.add_argument("--max-gap-size", type=int,   default=100000)
    return p.parse_args()


# BSJ detection ─────────────────────────────────────────────────────────────
def detect_bsj_paf(
    file_paf_gz,
    min_align_len=30,
    min_align_frac=0.35,
    max_bsj_gap=-3,
    max_bsj_ovlp=15,
    min_pident=0.95,
    min_m_quality=20,
    valid_coord_dist=10,
    max_gap_size=100000,
):
    PAF_DTYPES = {
        "q_name": "string",  # nullable string — read name
        "q_length": "int64",
        "q_start": "int64",
        "q_end": "int64",
        "s_strand": "string",  # only ever "+" or "-" — category is fastest for
        # the repeated equality checks (aln1["s_strand"] == "-")
        "s_chr_name": "string",  # small fixed set of chromosome/contig names,
        # repeated across millions of rows — category
        # avoids storing the same string over and over
        "s_length": "int64",
        "s_start": "int64",
        "s_end": "int64",
        "r_matches": "int64",
        "a_length": "int64",
        "m_quality": "int64",
        "x13": "string",  # optional SAM-style tags, e.g. "tp:A:P"
        "x14": "string",  # kept as nullable string, not numeric —
        "x15": "string",  # their format varies (tag:type:value) and
        "x16": "string",  # isn't meant to be parsed as a single number
        "x17": "string",
        "x18": "string",
    }

    paf_cols = [
        "q_name","q_length","q_start","q_end","s_strand","s_chr_name",
        "s_length","s_start","s_end","r_matches","a_length","m_quality",
        "x13","x14","x15","x16","x17","x18"
    ]
    paf_df = pd.read_csv(
        file_paf_gz,
        sep="\t",
        header=None,
        names=paf_cols,
        dtype=PAF_DTYPES, compression="gzip"
    )
    # paf_df = pd.read_csv(file_paf_gz, sep="\t", header=None, names=paf_cols)
    paf_df["pident"] = paf_df["r_matches"] / paf_df["a_length"]

    required_len = (paf_df["q_length"] * min_align_frac).clip(lower=min_align_len)

    paf_df = paf_df[
        (paf_df["a_length"] >= required_len)      &
        (paf_df["pident"]   >= min_pident)         &
        (paf_df["a_length"] <= paf_df["q_length"]) &
        (paf_df["a_length"] != paf_df["q_length"]) &
        (paf_df["m_quality"] >= min_m_quality)
    ].copy()

    logging.info(f"PAF rows after filter: {len(paf_df):,}")

    bsj_hits, val_hits = {}, {}

    for q_id, group in paf_df.groupby("q_name"):
        bsj_hits.setdefault(q_id, [])
        val_hits.setdefault(q_id, [])

        for s_id, sub_group in group.groupby("s_chr_name"):
            if len(sub_group) < 2:
                continue
            sub_group = sub_group.sort_values("q_start")

            for i in range(len(sub_group)):
                for j in range(i + 1, len(sub_group)):
                    aln1 = sub_group.iloc[i]
                    aln2 = sub_group.iloc[j]

                    if aln1["s_strand"] != aln2["s_strand"]:
                        continue
                    if abs(aln1["s_start"] - aln2["s_start"]) > max_gap_size:
                        continue

                    q_ovlp = aln1["q_end"] - aln2["q_start"]
                    if q_ovlp < max_bsj_gap or q_ovlp > max_bsj_ovlp:
                        continue

                    # ── FIXED: local variables, never mutate pandas rows ──────
                    ovl_len = int(q_ovlp) if q_ovlp > 0 else 0

                    if aln1["s_strand"] == "-":
                        if aln2["q_start"] > aln1["q_start"] and \
                           aln2["s_start"] > aln1["s_end"]:
                            s_end_corrected = int(aln2["s_end"]) - ovl_len
                            bsj_hits[q_id].append([
                                int(aln1["q_end"]), s_id,
                                int(aln1["s_start"]), s_end_corrected,
                                aln1["s_strand"], ovl_len
                            ])
                        else:
                            val_hits[q_id].append([int(aln1["q_end"]), s_id])
                    else:
                        if aln2["q_start"] > aln1["q_start"] and \
                           aln2["s_start"] < aln1["s_end"]:
                            s_start_corrected = int(aln2["s_start"]) + ovl_len
                            bsj_hits[q_id].append([
                                int(aln1["q_end"]), s_id,
                                s_start_corrected, int(aln1["s_end"]),
                                aln1["s_strand"], ovl_len
                            ])
                        else:
                            val_hits[q_id].append([int(aln1["q_end"]), s_id])

        for hit in bsj_hits[q_id]:
            is_bsj = all(
                abs(hit[0] - vh[0]) > valid_coord_dist
                for vh in val_hits[q_id]
            )
            if is_bsj:
                yield (
                    f"{q_id}\t{hit[0]}\t{hit[1]}\t"
                    f"{hit[2]}\t{hit[3]}\t{hit[4]}\t{hit[5]}\n"
                )


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


# ── Main ──────────────────────────────────────────────────────────────────────
def main():                               # ← FIXED: only ONE main()
    args = parse_args()

    Path(args.tsv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.fasta).parent.mkdir(parents=True, exist_ok=True)

    # Step 1 — write BSJ TSV
    n_bsj = 0
    with gzip.open(args.tsv, "wt") as out_tsv:
        out_tsv.write(
            "q_id\tbsj_coord\tchromosome\t"
            "s_start\ts_end\ts_strand\toverlap_trimmed_nt\n"
        )
        for line in detect_bsj_paf(
            args.paf,
                file_paf_gz,
                min_align_len=30,
                min_align_frac=0.35,
                max_bsj_gap=-3,
                max_bsj_ovlp=15,
                min_pident=0.95,
                min_m_quality=20,
                valid_coord_dist=10,
                max_gap_size=100000,
            ):
            out_tsv.write(line)
            n_bsj += 1

    logging.info(f"BSJ entries: {n_bsj:,} → {args.tsv}")
    if n_bsj == 0:
        logging.warning("No BSJ detected — check PAF and filter parameters")
        return

    bsj_coords = load_bsj_coordinates(args.bsj_tsv)

    # Step 2 — extract sequences
    extract_bsj_sequences(
        fastq_file     = args.fastq,
        bsj_coords     = bsj_coords,
        flank          = args.flank
    )
    logging.info("paf_parser.py complete ✓")


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":        # ← FIXED: top-level code moved here
    main()                        #   nothing runs on import