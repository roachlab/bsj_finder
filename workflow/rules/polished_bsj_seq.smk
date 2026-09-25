
# rule concat_trim_reads:
#     input:
#         os.path.join(out_dir,"trim_results","{sample}_R1_trim_rename.fastq.gz"),
#         os.path.join(out_dir,"trim_results","{sample}_R2_trim_rename.fastq.gz"),
#         os.path.join(out_dir, "trim_results", "{sample}_RS_trim_rename.fastq.gz"),
#     output:
#         os.path.join(out_dir,"trim_results","{sample}_trimmed.fastq.gz")
#     shell:
#         "cat {input} > {output}"

# Polishing uses ONLY the BSJ-supporting reads, pulled back out of the alignment BAM.
#
# bsj_cor.tsv.gz column 1 is the read id of every read that passed the chimeric
# geometry filters, and align_results/{sample}.bam still holds those records (it is
# already `-d SA` filtered, so it is a junction-read BAM). Selecting that id list and
# converting back to FASTQ gives the exact read set that supports the reference being
# polished -- on GBM sample 150 that is 133k reads instead of ~29.9M, a ~224x cut in
# alignment work, at no measurable cost in consensus accuracy (benchmarked against
# BSJ_test_dataset truth: mean edit distance to the 48 truth junction windows 3.92 vs
# 3.83 for the all-reads version, identical exact-match and <=3 counts).
#
# It also means this rule no longer consumes the temp() renamed reads. extract_bsj_seq
# is now their last consumer, so re-running polish (or anything downstream of it) no
# longer drags re-trimming and re-alignment back into the DAG.
rule polished_bsj_seq:
    input:
        ref_gz = os.path.join(out_dir, "dedup_results", "{sample}.bsj_dedup.fasta.gz"),
        bam    = os.path.join(out_dir, "align_results", "{sample}.bam"),
        cor    = os.path.join(out_dir, "bsj_results", "{sample}.bsj_cor.tsv.gz"),

    output:
        fasta = os.path.join(out_dir, "polished_results", "{sample}.polished.fasta"),
        fai   = os.path.join(out_dir, "polished_results", "{sample}.polished.fasta.fai")

    conda:
        os.path.join("..", "env", "polished.yaml")

    # Feeds minimap2 -t, samtools sort -@ and bcftools --threads below. Without this
    # the rule silently ran at threads=1.
    #
    # 4, not 8: measured on GBM sample 150 (68,083 refs, 133k BSJ reads) the whole
    # chain takes 20s / 14s / 13s at 1 / 4 / 8 threads. It plateaus after 4 because
    # the cost is bcftools mpileup walking the references serially, which does not
    # parallelise -- --threads only buys BGZF (de)compression. Spending the rest of
    # the budget on concurrent SAMPLES is strictly better: at --cores 40 this gives
    # 10 samples at once instead of 5.
    #
    # NB samtools sort -m is PER THREAD, so sort_mem (default 1G) x threads is the
    # per-sample ceiling -- raise sort_mem with this in mind.
    threads:
        4

    log:
        os.path.join("logs", "polished_bsj_seq_{sample}.log")

    params:
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

        # Pull just the BSJ-supporting reads out of the alignment BAM.
        #
        # --exclude-flags 0x900 drops secondary (0x100) and supplementary (0x800)
        # records. That is essential, not tidiness: a supplementary record is the
        # second junction arm and its SEQ is hard-clipped, so including it would feed
        # truncated sequence into the consensus. The primary record carries the full
        # read, and samtools fastq restores original orientation for reverse-strand
        # hits. -n keeps the read name as-is (no /1 /2 suffix); these were aligned
        # single-end, so everything is unpaired and lands on stdout.
        IDS="{params.workdir}/{wildcards.sample}.bsj_read_ids.txt"
        zcat {input.cor} | tail -n +2 | cut -f1 | sort -u > "$IDS"

        if [[ ! -s "$IDS" ]]; then
            # Reference is non-empty but nothing supports it -- leave it unpolished
            # rather than letting an empty pileup produce a bogus consensus.
            echo "[$(date)] No BSJ read ids for {wildcards.sample} — emitting reference unchanged" >> {log}
            cp "$REF" {output.fasta}
            samtools faidx {output.fasta} 2>> {log}
            rm -rf {params.workdir}
            exit 0
        fi

        samtools view -h -N "$IDS" --exclude-flags 0x900 {input.bam} 2>> {log} \
            | samtools fastq -n - 2>> {log} \
            | gzip -1 > "$READS"

        echo "[$(date)] BSJ read ids: $(wc -l < "$IDS") | reads extracted: $(zcat "$READS" | awk 'END{{print NR/4}}')" >> {log}
        echo "[$(date)] Aligning BSJ reads with minimap2" >> {log}
        
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
