# ═════════════════════════════════════════════════════════════════════════════
# Snakefile — circRNA BSJ Detection Pipeline
# ═════════════════════════════════════════════════════════════════════════════
import os
import re
from metasnek import fastq_finder

# Resolved against the Snakefile, NOT the working directory. A bare
# os.path.join("config", "config.yaml") is CWD-relative, so any run with
# `-d <dir>` outside the repo dies at parse time with "Workflow defines configfile
# config/config.yaml but it is not present or accessible" -- which is what forces
# every invocation to happen from the repo root. (READ_DIR below is normalised the
# same way.)
configfile: os.path.join(workflow.basedir, "config", "config.yaml")

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

# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE EXCLUSION
#
# Drop named samples without touching readDirectory. Excluding via a curated
# symlink directory (as the PRJEB58914 subsets do) would change every remaining
# sample's INPUT PATH, which is a Snakemake rerun trigger -- so a purely
# bookkeeping change would discard all banked alignment work. Naming them here
# leaves the paths identical and the existing outputs valid.
#
# The reason for each exclusion belongs in the config next to the name.
# ─────────────────────────────────────────────────────────────────────────────
EXCLUDE = set(config.get("exclude_samples") or [])
if EXCLUDE:
    unknown = EXCLUDE - set(SAMPLES)
    if unknown:
        # A typo here would silently exclude nothing, so fail instead.
        raise ValueError(
            f"exclude_samples names {sorted(unknown)}, which are not in {READ_DIR}. "
            f"Known samples: {SAMPLES}"
        )
    SAMPLES = [s for s in SAMPLES if s not in EXCLUDE]
    print(f"  excluded {len(EXCLUDE)}: {sorted(EXCLUDE)} -> {len(SAMPLES)} samples remain")
    if not SAMPLES:
        raise ValueError("exclude_samples removed every sample")


# ─────────────────────────────────────────────────────────────────────────────
# PAIRED-END / SINGLE-END PARTITION
#
# metasnek always emits an "R2" key but sets it to None for single-end samples,
# so truthiness -- not key presence -- is the test.
#
# Snakemake cannot vary the NUMBER of a rule's outputs per wildcard, so the trim
# and rename rules are split into _pe/_se variants: only then can R2/RS be real
# declared outputs rather than untracked params, which is what lets temp() reach
# them. Each variant is pinned to its own sample set via wildcard_constraints, so
# the two never compete for the same {sample}.
# ─────────────────────────────────────────────────────────────────────────────
PE_SAMPLES = [s for s in SAMPLES if samples[s].get("R2")]
SE_SAMPLES = [s for s in SAMPLES if not samples[s].get("R2")]
print(f"  paired-end: {len(PE_SAMPLES)} | single-end: {len(SE_SAMPLES)}")


def sample_regex(sample_list):
    """Regex matching exactly these sample names -- or nothing, if the list is empty.

    An empty join would be "", which matches ANY wildcard, so the _pe and _se rules
    would both claim every sample on a dataset that is entirely one or the other.
    "(?!)" is a never-satisfiable lookahead and is the correct empty case.
    """
    return "|".join(re.escape(s) for s in sample_list) if sample_list else r"(?!)"


def renamed_reads(wildcards):
    """Declared inputs for every rule that consumes renamed trim output.

    Unlike outputs, a rule's INPUTS may vary per wildcard via an input function, so
    the PE/SE difference is handled here instead of splitting each consumer. Rules
    use it as `unpack(renamed_reads)`; the shell bodies still reference the params
    copies of these paths, which keeps their existing `-s`/`-e` guards untouched.

    Every consumer of R2/RS must go through this, or temp() would free them as soon
    as the last rule Snakemake KNOWS about finishes -- while a later rule still
    needs them. Current consumers: minimap2_alig and extract_bsj_seq. (polished_bsj_seq
    dropped out when it switched to reading BSJ reads from the alignment BAM.)
    """
    base = os.path.join(out_dir, "trim_results", wildcards.sample)
    reads = {"r1": f"{base}_R1_trim_rename.fastq.gz"}
    if samples[wildcards.sample].get("R2"):
        reads["r2"] = f"{base}_R2_trim_rename.fastq.gz"
        reads["rs"] = f"{base}_RS_trim_rename.fastq.gz"
    return reads


# ─────────────────────────────────────────────────────────────────────────────
# MMSEQS2 BINARY SELECTION
#
# Default (mmseqs_static: false): the mmseqs rules use the conda env
# workflow/env/mmseqs2.yaml, like every other rule in the pipeline. That is the
# portable choice and is what a fresh checkout should do.
#
# mmseqs_static: true instead uses the pre-built SSE4.1 release binary vendored in
# tools/mmseqs/bin/ (from https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz) and
# attaches NO conda env to those rules. This is only needed on a host without AVX2:
# the conda/bioconda mmseqs2 package is AVX2-compiled and dies with SIGILL
# ("Illegal instruction") there. The runs in this working copy are on such a host
# (Xeon E5-2690 v2), so their configs set it -- see config/config_gbm.yaml and
# config/config_tassie_devil.yaml.
#
# The rules consume both values: MMSEQS_CONDA as their `conda:` (None = no env),
# and MMSEQS_BIN_DIR as a param prepended to PATH (empty = leave PATH alone).
# ─────────────────────────────────────────────────────────────────────────────
def as_bool(value):
    """Coerce a config value to bool, treating the usual YAML/CLI spellings as such.

    --config mmseqs_static=false arrives as the STRING "false", which is truthy --
    so a bare bool() here would silently turn the flag on when asked to turn it off.
    """
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


MMSEQS_STATIC = as_bool(config.get("mmseqs_static", False))
MMSEQS_BIN_DIR = os.path.join(workflow.basedir, "tools", "mmseqs", "bin") if MMSEQS_STATIC else ""
MMSEQS_CONDA = None if MMSEQS_STATIC else os.path.join(workflow.basedir, "workflow", "env", "mmseqs2.yaml")

if MMSEQS_STATIC:
    mmseqs_exe = os.path.join(MMSEQS_BIN_DIR, "mmseqs")
    if not os.access(mmseqs_exe, os.X_OK):
        # Fail here rather than let the rules fall through to whatever mmseqs is on
        # PATH -- on this host that is either nothing, or the AVX2 build that SIGILLs
        # hours into the run.
        raise FileNotFoundError(
            f"mmseqs_static is set but there is no executable mmseqs at {mmseqs_exe}. "
            "Fetch it with: curl -sL https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz "
            "| tar xz -C tools/  (note: tools/ is gitignored, so a fresh clone has "
            "neither the binary nor the tarball). Or set mmseqs_static: false to use "
            "the conda build."
        )
    print(f"mmseqs2: static SSE4.1 build ({mmseqs_exe})")
else:
    print(f"mmseqs2: conda env ({MMSEQS_CONDA})")


# ─────────────────────────────────────────────────────────────────────────────
# LIBRARY STRANDEDNESS
#
# Two separate things depend on this, and they are easy to conflate:
#
# 1. R2 IS ALWAYS REVERSE-COMPLEMENTED for paired samples (rename_read_ids_pe),
#    independent of this setting. Reads are aligned SINGLE-END, so an FR-oriented
#    fragment puts its mates on opposite genomic strands and the same junction gets
#    called twice, once per strand. Measured before the fix: 54.35% of GBM summary
#    rows and 5.60% of Tasmanian devil rows were +/- duplicate pairs, 98.9% of them
#    sharing overlap_trimmed_nt, with support split between the two entries. In 8
#    devil samples, 82.2% of both-strand coordinates had each mate on exactly one
#    strand -- the signature of this artefact. RC'ing R2 makes a fragment's mates
#    agree, collapsing the pair. Single-end samples are unaffected.
#
# 2. THIS SETTING decides how the resulting strand is INTERPRETED, and is applied
#    in sa_detect.py:
#      reverse    - R1 is antisense (dUTP / TruSeq Stranded). The junction strand is
#                   FLIPPED, because it otherwise follows R1's orientation.
#      forward    - R1 is sense (ligation-based). Reported as aligned.
#      unstranded - read orientation carries no strand information. Reported as
#                   aligned, but the strand column is NOT meaningful: R1 itself is
#                   randomly oriented, so RC'ing R2 cannot collapse everything and
#                   residual +/- duplicates are expected.
#
# Both real datasets measured so far are `reverse`, confirmed two independent ways
# on GBM sample 150: minimap2 ts:A splice-motif orientation (R1 96.2% antisense) and
# exonic reads vs GENCODE v44 gene strand (R1 97.2% antisense). The splice-motif
# method needs no annotation and is the one to use on genomes that lack one.
# ─────────────────────────────────────────────────────────────────────────────
STRANDEDNESS = str(config.get("library_strandedness", "auto")).strip().lower()
if STRANDEDNESS not in ("auto", "forward", "reverse", "unstranded"):
    raise ValueError(
        f"library_strandedness must be auto, forward, reverse or unstranded "
        f"(got {STRANDEDNESS!r})."
    )
print(f"library strandedness: {STRANDEDNESS}"
      + (" (detected per sample by detect_strandedness)" if STRANDEDNESS == "auto"
         else " (declared; detection still runs and fails on disagreement)"))


# ─────────────────────────────────────────────────────────────────────────────
# SEEDING PRESET
#
# minimap2 seed length (-k) and minimizer window (-w), chosen by name. Measured on
# BSJ_test_dataset (48 true BSJs; 8 threads, one job at a time):
#
#   sensitive  -k 11 -w 5   Dense seeding: recovers short and noisy junction arms. On the
#                           1.5 Mb test genome 44-48/48 BSJs in every library (100/150 bp,
#                           low/high error). On a 1 Gbp genome: 2.3 h for 36M read pairs.
#   fast       -k 21 -w 11  minimap2's `sr` values. 6-12x faster on a 1 Gbp genome (22 min
#                           for 36M pairs), all 48/48 on 150 bp low-error reads at any depth,
#                           but on 100 bp or high-error reads it loses low-coverage (1-3 read)
#                           circles (43-45/48) and reports 13-38% fewer junction reads; on
#                           small, short, noisy libraries strandedness auto-detection can see
#                           too few spliced reads, so declare library_strandedness there.
#
# Use `fast` for large genomes with >= 150 bp reads at decent depth; keep `sensitive` when
# reads are short or noisy, or when low-abundance circles matter. Explicit minimap2_k /
# minimap2_w in the config override the preset.
# ─────────────────────────────────────────────────────────────────────────────
SEED_PRESETS = {"sensitive": (11, 5), "fast": (21, 11)}
PRESET = str(config.get("preset", "sensitive")).strip().lower()
if PRESET not in SEED_PRESETS:
    raise ValueError(f"preset must be one of {sorted(SEED_PRESETS)} (got {PRESET!r}).")
_k, _w = SEED_PRESETS[PRESET]
config["minimap2_k"] = int(config.get("minimap2_k") or _k)
config["minimap2_w"] = int(config.get("minimap2_w") or _w)
print(f"seeding preset: {PRESET} -> minimap2 -k {config['minimap2_k']} -w {config['minimap2_w']}"
      + ("" if (config["minimap2_k"], config["minimap2_w"]) == (_k, _w) else " (minimap2_k/_w overridden)"))


include: "workflow/rules/trimnami_trim.smk"
include: "workflow/rules/strandedness.smk"
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

