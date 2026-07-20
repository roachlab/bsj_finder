#rule to run trimnami

rule trimnami_fastp:
    input:
        r1= lambda wc: samples[wc.sample]["R1"],
    output:
        r1= os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz"),
    params:
        r2_in = lambda wc: samples[wc.sample].get("R2", ""),
        r2_out = os.path.join(out_dir, "trim_results", "{sample}_R2_trim.fastq.gz"),
        orphan= os.path.join(out_dir, "trim_results", "{sample}_RS_trim.fastq.gz"),
    conda:
        os.path.join("..", "env", "trimnami.yaml")
    log:
        os.path.join("logs", "trimnami","{sample}.log")
    threads:
        8
    shell:
        """
         if [[ -n "{params.r2_in}" && -e "{params.r2_in}" ]]; then \\
         echo "[$(date)] Running PAIRED-END trimming for {wildcards.sample}" >> {log}
            fastp \\
                -i {input.r1} \\
                -I {params.r2_in} \\
                -o {output.r1} \\
                -O {params.r2_out} \\
                --unpaired1 {params.orphan} \\
                --unpaired2 {params.orphan} \\
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
        
        else \\
            echo "[$(date)] Running SINGLE-END trimming for {wildcards.sample}" >> {log}
            fastp \\
                -i {input.r1} \\
                -o {output.r1} \\
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
                2> {log}
        fi
        """



