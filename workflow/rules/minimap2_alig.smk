#rule to run trimnami

rule minimap2_alig:
    input:
        genome = config['genome'],
        r1 = os.path.join(out_dir, "trim_results", "{sample}_R1_trim_rename.fastq.gz")
    output:
        results = os.path.join(out_dir, "align_results", "{sample}.paf.gz")
    conda:
        os.path.join("..", "env", "minimap2.yaml")
    log:
        os.path.join("logs", "minimap2_alig.{sample}.log")
    threads:
        8
    params:
        r2_in = os.path.join(out_dir, "trim_results", "{sample}_R2_trim_rename.fastq.gz"),
        outdir = os.path.join(out_dir, "align_results", "{sample}.paf.gz")
    shell:
            """            
            if [[ -n "{params.r2_in}" && -e "{params.r2_in}" ]]; then
            minimap2 \\
                -x splice:sr \\
                --secondary=no \\
                -t {threads} \\
                {input.genome} \\
                {input.r1} \\
                {params.r2_in} \\
                2> {log} \\
                | gzip > {output.results}
            
            else
            
            minimap2 \\
                -x splice:sr \\
                --secondary=no \\
                -t {threads} \\
                {input.genome} \\
                {input.r1} \\
                2> {log} \\
                | gzip > {output.results}
            
            fi
            """

