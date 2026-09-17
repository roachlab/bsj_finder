# Align trimmed reads with minimap2 in SAM mode for chimeric (SA-tag) BSJ detection.
#
# All trimmed reads (R1, R2, orphans) are concatenated and aligned SINGLE-END with
# `-ax sr`, so every read is evaluated independently for the back-splice chimeric
# signature. minimap2 then reports each junction read's two arms as a primary
# alignment plus an SA:Z supplementary tag, which the detect step parses.
#
# The SAM stream is encoded to BAM on the fly and filtered to SA-bearing records only
# (`-d SA`): a read that maps end-to-end at one locus has no supplementary alignment and
# so no SA:Z tag, and sa_detect.py reconstructs both arms from the primary line's SA:Z
# string alone. Everything without an SA tag is therefore dead weight for detection --
# on real RNA-seq that is >99% of records. Supplementary records ARE the second junction
# arm, so they are kept for inspection in IGV.
#
# NOTE: this means align_results/{sample}.bam is a junction-read BAM, not a record of the
# whole library -- mapping rate, flagstat and linear/host-gene coverage cannot be computed
# from it. Drop `-d SA` if you need those.
#
# The BAM is left UNSORTED: detection is a single streaming pass over records, so a
# coordinate sort would be wasted work.

rule minimap2_alig:
    input:
        genome = config['genome'],
        r1 = os.path.join(out_dir, "trim_results", "{sample}_R1_trim_rename.fastq.gz")
    output:
        bam = os.path.join(out_dir, "align_results", "{sample}.bam")
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

        # Concatenate and align single-end; keep only SA-bearing (chimeric) records.
        zcat "${{READS[@]}}" \
            | minimap2 \
                -ax sr \
                -k 11 -w 5 \
                -t {threads} \
                {input.genome} \
                - \
                2> {log} \
            | samtools view -b -@ 2 -d SA -o {output.bam} -
        """
