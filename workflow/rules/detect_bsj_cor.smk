# Detect back-splice junctions from minimap2 SAM using the SA-tag / soft-clip parser.
# Consumes the gzipped SAM from minimap2_alig and emits one row per supporting read in
# the bsj_cor schema (q_id, bsj_coord, chromosome, s_start, s_end, s_strand,
# overlap_trimmed_nt) consumed by extract_bsj_seq. sa_detect.py is stdlib-only, so no
# conda environment is required.

rule detect_bsj_cor:
    input:
        sam = os.path.join(out_dir, "align_results", "{sample}.sam.gz")

    output:
        tsv = os.path.join(out_dir, "bsj_results", "{sample}.bsj_cor.tsv.gz")

    log:
        os.path.join("logs", "detect_bsj_cor.{sample}.log")

    params:
        script        = os.path.join(workflow.basedir, "workflow", "scripts", "sa_detect.py"),
        min_anchor    = config.get("min_align_len",     30),   # min aligned bp per arm
        min_mapq      = config.get("min_mapq",           1),
        max_gap_size  = config.get("max_gap_size",  100000),   # max genomic span of a circle
        max_read_gap  = config.get("sa_max_read_gap",   15),   # max unaligned gap between arms (read bp)
        max_read_ovlp = config.get("sa_max_read_ovlp",  15)    # max arm overlap / microhomology (read bp)

    shell:
        """
        mkdir -p "$(dirname {output.tsv})"
        zcat {input.sam} \
            | python3 {params.script} --sam - --out - \
                --min-anchor {params.min_anchor} \
                --min-mapq {params.min_mapq} \
                --max-gap-size {params.max_gap_size} \
                --max-read-gap {params.max_read_gap} \
                --max-read-ovlp {params.max_read_ovlp} \
                2> {log} \
            | gzip > {output.tsv}
        """
