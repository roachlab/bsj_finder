import os

rule gzip_tsv:
    input:
        "{path}.tsv"
    output:
        "{path}.tsv.gz"
    shell:
        "gzip {input}"


rule detect_bsj_cor:
    input:
        paf = os.path.join(out_dir, "align_results", "{sample}.paf.gz")

    output:
        tsv = os.path.join(out_dir, "bsj_results", "{sample}.bsj_cor.tsv")

    conda:
        os.path.join("..", "env", "biopython.yaml")
    log:
        os.path.join("logs", "bsj_detect", "{sample}.log")

    params:
        min_align_len = config.get("min_align_len", 50),
        max_bsj_gap = config.get("max_bsj_gap", -3),
        max_bsj_ovlp = config.get("max_bsj_ovlp", 15),
        min_pident = config.get("min_pident", 0.95),
        min_mapq = config.get("min_mapq", 20),
        valid_coord_dist = config.get("valid_coord_dist", 10),
        max_gap_size = config.get("max_gap_size", 100000)

    script:
        os.path.join("..", "scripts", "detect_bsj_cor.py")

