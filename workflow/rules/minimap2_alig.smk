# Align trimmed reads with minimap2 in SAM mode for chimeric (SA-tag) BSJ detection.
#
# All trimmed reads (R1, R2, orphans) are concatenated and aligned SINGLE-END with
# `-ax sr`, so every read is evaluated independently for the back-splice chimeric
# signature. minimap2 then reports each junction read's two arms as a primary
# alignment plus an SA:Z supplementary tag, which the detect step parses.

rule minimap2_alig:
    input:
        genome = config['genome'],
        r1 = os.path.join(out_dir, "trim_results", "{sample}_R1_trim_rename.fastq.gz")
    output:
        sam = os.path.join(out_dir, "align_results", "{sample}.sam.gz")
    conda:
        os.path.join("..", "env", "minimap2.yaml")
    log:
        os.path.join("logs", "minimap2_alig.{sample}.log")
    threads:
        8
    params:
        r2_in = os.path.join(out_dir, "trim_results", "{sample}_R2_trim_rename.fastq.gz"),
        orphan_in = os.path.join(out_dir, "trim_results", "{sample}_RS_trim_rename.fastq.gz")
    shell:
        """
        # Gather every trimmed read file that exists and has content.
        READS=("{input.r1}")
        if [ -s "{params.r2_in}" ]; then READS+=("{params.r2_in}"); fi
        if [ -s "{params.orphan_in}" ]; then READS+=("{params.orphan_in}"); fi

        # Concatenate and align single-end; emit gzipped SAM (carries SA:Z tags).
        zcat "${{READS[@]}}" \
            | minimap2 \
                -ax sr \
                -k 11 -w 5 \
                -t {threads} \
                {input.genome} \
                - \
                2> {log} \
            | gzip > {output.sam}
        """
