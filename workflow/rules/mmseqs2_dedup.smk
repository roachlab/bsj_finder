
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

    # mmseqs comes from the conda env by default; with `mmseqs_static: true` the
    # env is dropped (conda: None) and the vendored SSE4.1 binary is used instead --
    # see the MMSEQS2 BINARY SELECTION block in the Snakefile.
    conda: MMSEQS_CONDA

    # mmseqs defaults to every core on the box while snakemake reserves only 1, so
    # without this many concurrent samples each grab all 40. 8 matches the other rules.
    threads:
        8

    log:
        os.path.join("logs", "mmseqs2_dedup_{sample}.log")

    params:
        min_seq_id = config.get("mmseqs_min_seq_id", 0.95),
        min_cov = config.get("mmseqs_min_cov",    0.70),
        cov_mode = config.get("mmseqs_cov_mode",   0),
        workdir = os.path.join(out_dir, "dedup_results", "{sample}_mmseqs2_tmp"),# All intermediate files go into a per-sample working directory. Prevents file conflicts when multiple samples run in parallel
        mmseqs_bin_dir = MMSEQS_BIN_DIR,

    shell:
        """
        set -euo pipefail

        # The static SSE4.1 build, when enabled, takes precedence over any mmseqs
        # already on PATH (including the conda env's). Empty => leave PATH alone.
        if [[ -n "{params.mmseqs_bin_dir}" ]]; then
            export PATH="{params.mmseqs_bin_dir}:$PATH"
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
            --threads {threads} \
            2>> {log}

        mmseqs linclust \
            "$SEQDB" \
            "$CLUDB" \
            "$TMPDIR" \
            --min-seq-id {params.min_seq_id} \
            -c {params.min_cov}    \
            --cov-mode {params.cov_mode}   \
            --threads {threads} \
            2>> {log}
            
        mmseqs result2repseq \
            "$SEQDB" \
            "$CLUDB" \
            "$REPDB" \
            --threads {threads} \
            2>> {log}
            
        mmseqs convert2fasta \
            "$REPDB" \
            {output.fasta} \
            2>> {log}

        # mmseqs createtsv writes PLAIN TEXT, so it cannot write {output.tsv} directly
        # (declared .tsv.gz). Write it inside the workdir, then compress into place --
        # this matches the 0- and 1-sequence early-exit branches above, which already
        # gzip their output. Must happen before the workdir is removed below.
        mmseqs createtsv \
            "$SEQDB" \
            "$SEQDB" \
            "$CLUDB" \
            "{params.workdir}/cluster.tsv" \
            --threads {threads} \
            2>> {log}

        gzip -c "{params.workdir}/cluster.tsv" > {output.tsv}
                  
        rm -rf {params.workdir}

        # Summary 
        OUTPUT_COUNT=$( grep -c ">" {output.fasta} || echo 0)
        # REMOVED=$((SEQ_COUNT - OUTPUT_COUNT))
        # echo "Input sequences  : $SEQ_COUNT"  2>> {log}
        echo "Unique sequences : $OUTPUT_COUNT" 2>> {log}
        # echo "Removed          : $REMOVED"      2>> {log}
        
        """


