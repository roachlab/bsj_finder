
# rule concat_trim_reads:
#     input:
#         os.path.join(out_dir,"trim_results","{sample}_R1_trim_rename.fastq.gz"),
#         os.path.join(out_dir,"trim_results","{sample}_R2_trim_rename.fastq.gz"),
#         os.path.join(out_dir, "trim_results", "{sample}_RS_trim_rename.fastq.gz"),
#     output:
#         os.path.join(out_dir,"trim_results","{sample}_trimmed.fastq.gz")
#     shell:
#         "cat {input} > {output}"

# todo: maybe update polish step to use only BSJ reads
rule polished_bsj_seq:
    input:
        ref_gz = os.path.join(out_dir, "dedup_results", "{sample}.bsj_dedup.fasta.gz"),
        r1  = os.path.join(out_dir, "trim_results",  "{sample}_R1_trim_rename.fastq.gz")

    output:
        fasta = os.path.join(out_dir, "polished_results", "{sample}.polished.fasta"),
        fai   = os.path.join(out_dir, "polished_results", "{sample}.polished.fasta.fai")

    conda:
        os.path.join("..", "env", "polished.yaml")

    log:
        os.path.join("logs", "polished_bsj_seq_{sample}.log")

    params:
        r2_in = os.path.join(out_dir, "trim_results", "{sample}_R2_trim_rename.fastq.gz"),
        orphan = os.path.join(out_dir, "trim_results", "{sample}_RS_trim_rename.fastq.gz"),
        outdir = os.path.join(out_dir, "polished_results", "{sample}.polished.fasta"),
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
        READS="{params.workdir}/{wildcards.sample}.reads.fastq.gz"
        BAM="{params.workdir}/{wildcards.sample}.sorted.bam"
        RAW_VCF="{params.workdir}/{wildcards.sample}.raw.vcf.gz"
        NORM_VCF="{params.workdir}/{wildcards.sample}.norm.vcf.gz"
        FILT_VCF="{params.workdir}/{wildcards.sample}.filtered.vcf.gz"
        POLISHED="{params.workdir}/{wildcards.sample}.polished.fasta"

        gzip -c -d {input.ref_gz} > "$REF" #Decompress reference
        
        #Guard: check reference has real content before running anything
        
        REF_SEQ_COUNT=$(grep -c "^>" "$REF" || echo 0)
        REF_SEQ_COUNT=$(echo "$REF_SEQ_COUNT" | tr -d '[:space:]')
        echo "[$(date)] Reference sequence count: $REF_SEQ_COUNT" >> {log}
        
        if [[ "$REF_SEQ_COUNT" -eq 0 ]]; then
            echo "[$(date)] No BSJ reference sequences for {wildcards.sample} — writing empty outputs" >> {log}
            touch {output.fasta}
            touch {output.fai}
            rm -rf {params.workdir}
            exit 0
        fi

        #Align raw reads to deduplicated BSJ sequences     
        #Build combined reads file — only include sources that actually exist and contain data, so single-end samples work automatically.
        # Build an array — each file is a separate, distinct element
        
        READ_FILES=("{input.r1}")
        if [[ -s "{params.r2_in}" ]]; then
            READ_FILES+=("{params.r2_in}")
        fi
        if [[ -s "{params.orphan}" ]]; then
            READ_FILES+=("{params.orphan}")
        fi
        cat "${{READ_FILES[@]}}" > "$READS"

        
        echo "[$(date)] Aligning reads with minimap2: ${{READ_FILES[@]}}" >> {log}
        
        minimap2 -x sr -a -t {threads} "$REF" "$READS" 2>> {log} \
            | samtools view -h -q {params.min_mapq} -m {params.min_match} \
            | samtools sort -@ {threads} -m {params.sort_mem} \
            > "$BAM" 2>> {log}
        

        samtools index "$BAM" 2>> {log} 

        samtools faidx "$REF" 2>> {log} 
   
        bcftools mpileup -Ou -f "$REF" --threads {threads} "$BAM" 2>> {log} \
            | bcftools call -mv -Oz --threads {threads} -o "$RAW_VCF" 2>> {log}

        bcftools index "$RAW_VCF" 2>> {log}

        bcftools norm -f "$REF" -Oz -o "$NORM_VCF" "$RAW_VCF" 2>> {log}

        bcftools index "$NORM_VCF" 2>> {log}

        bcftools filter -Oz -o "$FILT_VCF" \
            -e 'QUAL<{params.min_qual} || INFO/DP<{params.min_dp}' \
            "$NORM_VCF" 2>> {log}
            
        bcftools index "$FILT_VCF" 2>> {log}

        echo "[$(date)] Generating polished consensus..." >> {log}
        bcftools consensus -f "$REF" "$FILT_VCF" > "$POLISHED" 2>> {log}

        mv "$POLISHED" {output.fasta}
        samtools faidx {output.fasta} 2>> {log}

        # Summary
        INPUT_SEQS=$(grep -c "^>" "$REF"      || true)
        OUTPUT_SEQS=$(grep -c "^>" "{output.fasta}" || true)
        RAW_VARIANTS=$(bcftools stats "$RAW_VCF"  | grep "number of records:" | head -1 | cut -f4)
        FILT_VARIANTS=$(bcftools stats "$FILT_VCF" | grep "number of records:" | head -1 | cut -f4)

        echo "=================================" >> {log}
        echo "Reference sequences : $INPUT_SEQS"    >> {log}
        echo "Polished sequences  : $OUTPUT_SEQS"   >> {log}
        echo "Raw variants        : $RAW_VARIANTS"  >> {log}
        echo "Filtered variants   : $FILT_VARIANTS" >> {log}
        echo "=================================" >> {log}

        rm -rf {params.workdir} #Clean up per-sample working directory
        """
