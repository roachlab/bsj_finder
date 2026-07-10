import os

rule bsj_cor_fasta:
    input:
        tsv   = os.path.join(out_dir, "bsj_results",  "{sample}.bsj_cor.tsv.gz"),
        fastq = os.path.join(out_dir, "trim_results",  "{sample}_{R12}_trim.fastq.gz")
    output:
        fasta = os.path.join(out_dir, "bsj_results", "{sample}_{R12}.bsj_sequences.fasta.gz")
    wildcard_constraints:
        R12=r'R[12]'
    conda:
        os.path.join("..", "env", "biopython.yaml")
    log:
        os.path.join("logs", "bsj_extract", "{sample}_{R12}.log")
    params:
        flank       = config.get("bsj_flank",      50),   # bases each side of BSJ
        read_suffix = config.get("read_suffix",    "/1")  # suffix on FASTQ read IDs
    script:
        os.path.join("..", "scripts", "extract_bsj_seq.py")
