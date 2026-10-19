# Detect back-splice junctions from the minimap2 BAM using the SA-tag / soft-clip parser.
# Consumes the unsorted, SA-filtered BAM from minimap2_alig and emits one row per
# supporting read in the bsj_cor schema (q_id, bsj_coord, chromosome, s_start, s_end,
# s_strand, overlap_trimmed_nt, mapq) consumed by extract_bsj_seq. sa_detect.py is stdlib-only
# and reads SAM text on stdin, so samtools streams the BAM back out as SAM for it; the
# shared minimap2 env supplies samtools.
#
# The `-d SA` filter already happened upstream, so the parser only ever sees chimeric
# records. `-F 0x900` additionally drops secondary/supplementary lines HERE rather than
# on disk: sa_detect.py anchors only on the primary line and rebuilds the second arm from
# that line's SA:Z string, so supplementary records are ~50% of the BAM and 100% wasted
# Python parsing -- but they are also the second junction arm, so they stay in the BAM
# for IGV. sa_detect.py keeps its own flag/SA checks so it remains correct standalone.
#
# This was a `samtools | python3 | gzip` shell pipeline; the samtools and gzip ends now
# live inside the script so `script:` can be used. That matters because Snakemake hashes
# a script: file as part of the rule's code, whereas the old `python3 {params.script}`
# form made it an untracked param string -- edits to the detector never triggered a
# rerun. Streaming is unchanged (samtools is still a subprocess read lazily).

rule detect_bsj_cor:
    input:
        bam = os.path.join(out_dir, "align_results", "{sample}.bam"),
        # An INPUT, not a param: params are resolved during DAG construction, before
        # detect_strandedness has run, so they cannot carry a computed verdict.
        strandedness = os.path.join(out_dir, "qc_results", "{sample}.strandedness.tsv")

    output:
        tsv = os.path.join(out_dir, "bsj_results", "{sample}.bsj_cor.tsv.gz")

    conda:
        os.path.join("..", "env", "minimap2.yaml")

    log:
        os.path.join("logs", "detect_bsj_cor.{sample}.log")

    threads:
        2

    params:
        min_anchor    = config.get("min_align_len",     30),   # min aligned bp per arm
        min_mapq      = config.get("min_mapq",           1),
        max_gap_size  = config.get("max_gap_size",  100000),   # max genomic span of a circle
        max_read_gap  = config.get("sa_max_read_gap",   15),   # max unaligned gap between arms (read bp)
        max_read_ovlp = config.get("sa_max_read_ovlp",  15)    # max arm overlap / microhomology (read bp)

    script:
        os.path.join("..", "scripts", "sa_detect.py")
