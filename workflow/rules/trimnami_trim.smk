#rule to run trimnami

rule trimnami_trim:
    input:
        r1 = lambda wc: samples[wc.sample]["R1"],
        r2 = lambda wc: samples[wc.sample]["R2"]

    output:
        r1 = os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz"),
        r2= os.path.join(out_dir, "trim_results", "{sample}_R2_trim.fastq.gz")
    conda:
        os.path.join("..", "env", "trimnami.yaml")

    log:
        os.path.join("logs", "trimnami","{sample}.log")

    params:
        outdir = os.path.join(out_dir, "trim_results")

    shell:
        """
        trimnami run \
            --reads {input} \
            --output {params.outdir} \
            2> {log}
        """

