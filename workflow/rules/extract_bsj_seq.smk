import os

rule bsj_cor_fasta:
    input:
        # see renamed_reads() in the Snakefile -- adds r2/rs for paired samples only
        unpack(renamed_reads),
        tsv = os.path.join(out_dir, "bsj_results",  "{sample}.bsj_cor.tsv.gz"),
    output:
        fasta = os.path.join(out_dir, "bsj_results", "{sample}_bsj_sequences.fasta.gz")
    conda:
        os.path.join("..", "env", "biopython.yaml")
    log:
        os.path.join("logs", "bsj_cor_fasta_{sample}.log")
    params:
        r2_in = os.path.join(out_dir, "trim_results", "{sample}_R2_trim_rename.fastq.gz"),
        orphan_in = os.path.join(out_dir, "trim_results", "{sample}_RS_trim_rename.fastq.gz"),
        fasta = os.path.join(out_dir, "bsj_results", "{sample}_bsj_sequences.fasta.gz"),
        flank       = config.get("bsj_flank",      50),   # bases each side of BSJ
        #read_suffix = config.get("read_suffix",    "/1")  # suffix on FASTQ read IDs
    script:
        os.path.join("..", "scripts", "extract_bsj_seq.py")
