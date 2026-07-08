import os

rule concat_trim_reads:
    input:
        os.path.join(out_dir,"trim_results","{sample}_R1_trim.fastq.gz"),
        os.path.join(out_dir,"trim_results","{sample}_R2_trim.fastq.gz"),
    output:
        os.path.join(out_dir,"trim_results","{sample}_trim.fastq.gz")
    shell:
        "cat {input} > {output}"


rule polished_bsj_seq:
    input:
        ref_gz = os.path.join(out_dir, "dedup_results", "{sample}.bsj_dedup.fasta.gz"),
        reads  = os.path.join(out_dir, "trim_results",  "{sample}_trim.fastq.gz")

    output:
        fasta = os.path.join(out_dir, "polished_results", "{sample}.polished.fasta"),
        fai   = os.path.join(out_dir, "polished_results", "{sample}.polished.fasta.fai")

    conda:
        os.path.join("..", "env", "polished.yaml")

    log:
        os.path.join("logs", "polish", "{sample}.log")

    params:
        min_mapq   = config.get("polish_min_mapq",   20),   # samtools view -q
        min_match  = config.get("polish_min_match",  75),   # samtools view -m
        min_qual   = config.get("bcftools_min_qual", 20),   # variant QUAL filter
        min_dp     = config.get("bcftools_min_dp",   3),    # variant DP filter
        sort_mem   = config.get("sort_mem", "1G"), # samtools sort -m
        workdir    = os.path.join(out_dir, "polished_results", "{sample}_polish_tmp") # Per-sample working directory — avoids filename clashes across samples

    shell:
        """
        set -euo pipefail

        mkdir -p {params.workdir}

        REF="{params.workdir}/{wildcards.sample}.ref.fasta"
        BAM="{params.workdir}/{wildcards.sample}.sorted.bam"
        RAW_VCF="{params.workdir}/{wildcards.sample}.raw.vcf.gz"
        NORM_VCF="{params.workdir}/{wildcards.sample}.norm.vcf.gz"
        FILT_VCF="{params.workdir}/{wildcards.sample}.filtered.vcf.gz"
        POLISHED="{params.workdir}/{wildcards.sample}.polished.fasta"

        # Step 1: Decompress reference ──────────────────────────────────────
        echo "[$(date)] Step 1: Decompressing reference..." >> {log}
        zcat {input.ref_gz} > "$REF"

        # Step 2: Align raw reads to deduplicated BSJ sequences ────────────
        echo "[$(date)] Step 2: Aligning reads with minimap2..." >> {log}
        minimap2 -x sr -a -t {threads} "$REF" {input.reads} 2>> {log} \
            | samtools view -h -q {params.min_mapq} -m {params.min_match} \
            | samtools sort -@ {threads} -m {params.sort_mem} \
            > "$BAM" 2>> {log}

        # Step 3: Index BAM file ─────────────────────────────────────────────
        echo "[$(date)] Step 3: Indexing BAM..." >> {log}
        samtools index "$BAM" 2>> {log}

        # Step 4: Index reference FASTA ─────────────────────────────────────
        echo "[$(date)] Step 4:Indexing reference FASTA..." >> {log}
        samtools faidx "$REF" 2>> {log}

        # Step 5: Call variants ──────────────────────────────────────────────
        echo "[$(date)] Step 5:Calling variants..." >> {log}
        bcftools mpileup -Ou -f "$REF" --threads {threads} "$BAM" 2>> {log} \
            | bcftools call -mv -Oz --threads {threads} -o "$RAW_VCF" 2>> {log}

        # Step 6: Index raw VCF ──────────────────────────────────────────────
        echo "[$(date)] Step 6: Indexing raw VCF..." >> {log}
        bcftools index "$RAW_VCF" 2>> {log}

        # Step 7: Normalise variants ──────────────────────────────────────────
        echo "[$(date)] Step 7: Normalising variants..." >> {log}
        bcftools norm -f "$REF" -Oz -o "$NORM_VCF" "$RAW_VCF" 2>> {log}

        # Step 8: Index normalised VCF ───────────────────────────────────────
        echo "[$(date)] Step 8: Indexing normalised VCF..." >> {log}
        bcftools index "$NORM_VCF" 2>> {log}

        # Step 9: Filter low-quality variants ────────────────────────────────
        echo "[$(date)] Step 9: Filtering variants (QUAL>={params.min_qual}, DP>={params.min_dp})..." >> {log}
        bcftools filter -Oz -o "$FILT_VCF" \
            -e 'QUAL<{params.min_qual} || INFO/DP<{params.min_dp}' \
            "$NORM_VCF" 2>> {log}

        # Step 10: Index filtered VCF ────────────────────────────────────────
        echo "[$(date)] Step 10: Indexing filtered VCF..." >> {log}
        bcftools index "$FILT_VCF" 2>> {log}

        # Step 11: Generate polished consensus FASTA ─────────────────────────
        echo "[$(date)] Step 11: Generating polished consensus..." >> {log}
        bcftools consensus -f "$REF" "$FILT_VCF" > "$POLISHED" 2>> {log}

        # ── Compress and index final polished FASTA ────────────────────────────
        mv "$POLISHED" > {output.fasta}
        samtools faidx {output.fasta} 2>> {log}

        # ── Summary ─────────────────────────────────────────────────────────────
        INPUT_SEQS=$(grep -c "^>" "$REF"      || true)
        OUTPUT_SEQS=$(grep -c "^>" "$POLISHED" || true)
        RAW_VARIANTS=$(bcftools stats "$RAW_VCF"  | grep "number of records:" | head -1 | cut -f4)
        FILT_VARIANTS=$(bcftools stats "$FILT_VCF" | grep "number of records:" | head -1 | cut -f4)

        echo "=================================" >> {log}
        echo "Reference sequences : $INPUT_SEQS"    >> {log}
        echo "Polished sequences  : $OUTPUT_SEQS"   >> {log}
        echo "Raw variants        : $RAW_VARIANTS"  >> {log}
        echo "Filtered variants   : $FILT_VARIANTS" >> {log}
        echo "=================================" >> {log}

        # ── Clean up per-sample working directory ──────────────────────────────
        rm -rf {params.workdir}
        """
