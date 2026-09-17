# ═════════════════════════════════════════════════════════════════════════════
# Snakefile — circRNA BSJ Detection Pipeline
# ═════════════════════════════════════════════════════════════════════════════
import os
from metasnek import fastq_finder

configfile: os.path.join("config", "config.yaml")

out_dir = config.get("outdir", "results")

READ_DIR = config.get("readDirectory")
if not READ_DIR:
    raise ValueError("Missing 'readDirectory' in config.yaml")

if not os.path.isabs(READ_DIR):
    READ_DIR = os.path.join(workflow.basedir, READ_DIR)

if not os.path.isdir(READ_DIR):
    raise NotADirectoryError(f"Sample directory not found: {READ_DIR}")

samples = fastq_finder.parse_samples_to_dictionary(READ_DIR)

if not samples:
    raise ValueError(f"No FASTQ samples found in: {READ_DIR}")

SAMPLES = sorted(samples.keys())
print(f"Samples detected ({len(SAMPLES)}): {SAMPLES}")


include: "workflow/rules/trimnami_trim.smk"
include: "workflow/rules/minimap2_alig.smk"
include: "workflow/rules/detect_bsj_cor.smk"
include: "workflow/rules/extract_bsj_seq.smk"
include: "workflow/rules/mmseqs2_dedup.smk"
include: "workflow/rules/polished_bsj_seq.smk"
include: "workflow/rules/pool_dedup.smk"


# ─────────────────────────────────────────────────────────────────────────────
# FINAL TARGETS
# ─────────────────────────────────────────────────────────────────────────────
rule all:
    input:
        expand(
            os.path.join(out_dir, "polished_results", "{sample}.polished.fasta"),
            sample=SAMPLES
        ),
        expand(
            os.path.join(out_dir, "polished_results", "{sample}.polished.fasta.fai"),
            sample=SAMPLES
        ),
        expand(
            os.path.join(out_dir, "dedup_results", "{sample}.bsj_cluster.tsv.gz"),
            sample=SAMPLES
        ),
        # Cross-sample pools: raw BSJ pool + final polished BSJ pool
        os.path.join(out_dir, "pool_results", "all_samples.bsj_pool_dedup.fasta"),
        os.path.join(out_dir, "pool_results", "all_samples.bsj_pool_cluster.tsv"),
        os.path.join(out_dir, "pool_results", "all_samples.polished_pool_dedup.fasta"),
        os.path.join(out_dir, "pool_results", "all_samples.polished_pool_cluster.tsv"),
        os.path.join(out_dir, "pool_results", "all_samples.bsj_summary.tsv")

