
rule gzip_fasta:
    input:
        "{path}.fasta"
    output:
        "{path}.fasta.gz"
    shell:
        "gzip -1 {input}"


rule mmseqs2_dedup:
    input:
        fasta = os.path.join(out_dir, "bsj_results", "{sample}_bsj_sequences.fasta.gz")

    output:
        fasta = os.path.join(out_dir, "dedup_results", "{sample}.bsj_dedup.fasta"),
        tsv = os.path.join(out_dir, "dedup_results", "{sample}.bsj_cluster.tsv.gz")

    # NOTE: no conda env here on this host. The conda/bioconda mmseqs2 build is
    # AVX2-compiled, but this machine (Xeon E5-2690 v2) has no AVX2 and the
    # binary crashes with SIGILL ("Illegal instruction"). Instead we use the
    # static SSE4.1 release binary shipped in tools/mmseqs/ (downloaded from
    # https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz), which runs on this CPU.
    # The static binary is self-contained; only standard coreutils are needed.

    log:
        os.path.join("logs", "mmseqs2_dedup_{sample}.log")

    params:
        min_seq_id = config.get("mmseqs_min_seq_id", 0.95),
        min_cov = config.get("mmseqs_min_cov",    0.70),
        cov_mode = config.get("mmseqs_cov_mode",   0),
        workdir = os.path.join(out_dir, "dedup_results", "{sample}_mmseqs2_tmp"),# All intermediate files go into a per-sample working directory. Prevents file conflicts when multiple samples run in parallel

    shell:
        """
        set -euo pipefail

        # Prefer the SSE4.1 static mmseqs (see note above); fall back to any
        # mmseqs already on PATH if the tools/ dir is missing.
        if [[ -x "tools/mmseqs/bin/mmseqs" ]]; then
            export PATH="$(pwd)/tools/mmseqs/bin:$PATH"
        fi

        # Create per-sample working directory
        mkdir -p {params.workdir}
        TMPDIR="{params.workdir}/tmp"
        mkdir -p "$TMPDIR"
        SEQDB="{params.workdir}/seqDB"
        CLUDB="{params.workdir}/cluDB"
        REPDB="{params.workdir}/repDB"
        
        echo "[$(date)] START dedup: {wildcards.sample}" > {log}
        
        SEQ_COUNT=$(gunzip -c {input.fasta} | grep -c "^>" || echo 0)
        SEQ_COUNT=$(echo "$SEQ_COUNT" | tr -d '[:space:]')
        echo "[$(date)] Input sequence count: $SEQ_COUNT" >> {log}

        if [[ "$SEQ_COUNT" -eq 0 ]]; then
            echo "[$(date)] No sequences found — writing empty outputs" >> {log}
            gunzip -c {input.fasta} > {output.fasta}
            printf "query\\ttarget\\n" | gzip > {output.tsv}
            rm -rf {params.workdir}
            exit 0
        fi

        if [[ "$SEQ_COUNT" -eq 1 ]]; then
            echo "[$(date)] Only 1 sequence — skipping clustering" >> {log}
            gunzip -c {input.fasta} > {output.fasta}
            gunzip -c {input.fasta} | grep "^>" | sed 's/^>//' | \
                awk '{{print $1"\\t"$1}}' | gzip > {output.tsv}
            rm -rf {params.workdir}
            exit 0
        fi
                
        mmseqs createdb \
            {input.fasta} \
            "$SEQDB" \
            2>> {log}

        mmseqs linclust \
            "$SEQDB" \
            "$CLUDB" \
            "$TMPDIR" \
            --min-seq-id {params.min_seq_id} \
            -c {params.min_cov}    \
            --cov-mode {params.cov_mode}   \
            2>> {log}
            
        mmseqs result2repseq \
            "$SEQDB" \
            "$CLUDB" \
            "$REPDB" \
            2>> {log}
            
        mmseqs convert2fasta \
            "$REPDB" \
            {output.fasta} \
            2>> {log}

        mmseqs createtsv \
            "$SEQDB" \
            "$SEQDB" \
            "$CLUDB" \
            {output.tsv} \
            2>> {log}
                  
        rm -rf {params.workdir}

        # Summary 
        OUTPUT_COUNT=$( grep -c ">" {output.fasta} || echo 0)
        # REMOVED=$((SEQ_COUNT - OUTPUT_COUNT))
        # echo "Input sequences  : $SEQ_COUNT"  2>> {log}
        echo "Unique sequences : $OUTPUT_COUNT" 2>> {log}
        # echo "Removed          : $REMOVED"      2>> {log}
        
        """


