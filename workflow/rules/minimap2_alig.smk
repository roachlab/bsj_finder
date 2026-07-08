#rule to run trimnami
import os
rule minimap2_alig:
    input:
        genome = config['genome'],
        r1 = os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz"),
        r2 = os.path.join(out_dir, "trim_results", "{sample}_R2_trim.fastq.gz")
    output:
        results = os.path.join(out_dir, "align_results", "{sample}.paf.gz")
    conda:
        os.path.join("..", "env", "minimap2.yaml")
    log:
        os.path.join("logs", "trimnami", "{sample}.log")
    params:
         outdir = os.path.join(out_dir, "align_results", "{sample}.paf.gz")
    shell:
            """
            minimap2 \
                -x sr \
                -t {threads} \
                {input.genome} \
                {input.r1} \
                {input.r2} \
                2> {log} \
            | gzip > {params.outdir}
            """