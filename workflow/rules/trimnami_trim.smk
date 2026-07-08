#rule to run trimnami

rule trimnami_fastp_paired_end:
    """Read trimming with fastp for paired reads"""
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
    threads:
        8
    shell:
        """
        fastp \\
            -i {input.r1} \\
            -I {input.r2} \\
            -o {output.r1} \\
            -O {output.r2} \\
            -z 1 \\
            --thread {threads} \\
            --qualified_quality_phred 15 \\
            --length_required 90 \\
            --cut_tail \\
            --cut_tail_window_size 25 \\
            --cut_tail_mean_quality 15 \\
            --dedup \\
            --dup_calc_accuracy 4 \\
            --trim_poly_x \\
            --detect_adapter_for_pe \\
            2> {log}
        """

# rule trimnami_trim:
#     input:
#         r1 = lambda wc: samples[wc.sample]["R1"],
#         r2 = lambda wc: samples[wc.sample]["R2"]
#
#     output:
#         r1 = os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz"),
#         r2= os.path.join(out_dir, "trim_results", "{sample}_R2_trim.fastq.gz")
#
#     conda:
#         os.path.join("..", "env", "trimnami.yaml")
#
#     log:
#         os.path.join("logs", "trimnami","{sample}.log")
#
#     params:
#         outdir = os.path.join(out_dir, "trim_results")
#
#     shell:
#         """
#         trimnami run \
#             --reads {input} \
#             --output {params.outdir} \
#             2> {log}
#         """

