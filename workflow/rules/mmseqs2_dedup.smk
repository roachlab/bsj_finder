import os

rule concat_bsj_extract:
    input:
        os.path.join(out_dir, "bsj_results", "{sample}_R1.bsj_sequences.fasta.gz"),
        os.path.join(out_dir, "bsj_results", "{sample}_R2.bsj_sequences.fasta.gz")
    output:
        os.path.join(out_dir,"bsj_results","{sample}_concat.bsj_sequences.fasta.gz")
    shell:
        "cat {input} > {output}"



rule mmseqs2_dedup:
    input:
        fasta = os.path.join(out_dir, "bsj_results", "{sample}_concat.bsj_sequences.fasta.gz")

    output:
        fasta = os.path.join(out_dir, "dedup_results", "{sample}.bsj_dedup.fasta.gz"),
        tsv = os.path.join(out_dir, "dedup_results", "{sample}.bsj_cluster.tsv.gz")

    conda:
        os.path.join("..", "env", "mmseqs2.yaml")

    log:
        os.path.join("logs", "dedup", "{sample}.log")

    params:
        min_seq_id = config.get("mmseqs_min_seq_id", 0.95),
        min_cov = config.get("mmseqs_min_cov",    0.70),
        cov_mode = config.get("mmseqs_cov_mode",   0),
        # All intermediate files go into a per-sample working directory. Prevents file conflicts when multiple samples run in parallel
        workdir = os.path.join(out_dir, "dedup_results", "{sample}_mmseqs2_tmp"),

    shell:
        """
        set -euo pipefail

        # Create per-sample working directory ───────────────────────────────
        mkdir -p {params.workdir}
        WORKDIR="{params.workdir}"
        TMPDIR="{params.workdir}/tmp"
        mkdir -p "$TMPDIR"

        # Intermediate database names (all inside workdir) ─────────────────
        SEQDB="{params.workdir}/seqDB"
        CLUDB="{params.workdir}/cluDB"
        REPDB="{params.workdir}/repDB"
        
        echo "[$(date)] START dedup: {wildcards.sample}" > {log}
        
        # Step 1. Decompress input FASTA ────────────────────────────────────────────
        echo "[$(date)] Step 1: Decompressing input FASTA..." >> {log}
        gunzip -c {input.fasta} > "$WORKDIR/{wildcards.sample}_input.fasta" 2>> {log}

        # Step 2: Create MMseqs2 database from input FASTA ─────────────────
        echo "[$(date)] Step 2: Creating sequence database....." 2>> {log}
        mmseqs createdb \
            {input.fasta} \
            "$SEQDB" \
            2>> {log}

        # Step 3: Cluster with linclust ─────────────────────────────────────
        echo "[$(date)] Step 3: Clustering sequences..." 2>> {log}
        mmseqs linclust \
            "$SEQDB" \
            "$CLUDB" \
            "$TMPDIR" \
            --min-seq-id {params.min_seq_id} \
            -c {params.min_cov}    \
            --cov-mode {params.cov_mode}   \
            2>> {log}

        # Step 4: Extract one representative per cluster ────────────────────
        echo "[$(date)] Step 4: Extracting representatives..." 2>> {log}
        mmseqs result2repseq \
            "$SEQDB" \
            "$CLUDB" \
            "$REPDB" \
            2>> {log}

        # Step 5: Convert representative database → FASTA ──────────────────
        echo "[$(date)] Step 5: Converting to FASTA..." 2>> {log}
        mmseqs convert2fasta \
            "$REPDB" \
            {output.fasta} \
            2>> {log}

        # Step 6: Export cluster membership table ───────────────────────────
        echo "[$(date)] Step 6: Writing cluster membership table..." 2>> {log}
        mmseqs createtsv \
            "$SEQDB" \
            "$SEQDB" \
            "$CLUDB" \
            {output.tsv} \
            2>> {log}

        # Clean up intermediate files ───────────────────────────────────────
        rm -rf {params.workdir}

        # Summary ───────────────────────────────────────────────────────────
        INPUT_COUNT=$(grep -c ">" {input.fasta}  || true)
        OUTPUT_COUNT=$(grep -c ">" {output.fasta} || true)
        REMOVED=$((INPUT_COUNT - OUTPUT_COUNT))
        echo "Input sequences  : $INPUT_COUNT"  2>> {log}
        echo "Unique sequences : $OUTPUT_COUNT" 2>> {log}
        echo "Removed          : $REMOVED"      2>> {log}

        echo "Done. Output: {output.fasta}" 2>> {log}
        """


