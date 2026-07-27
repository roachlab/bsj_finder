
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

    conda:
        os.path.join("..", "env", "mmseqs2.yaml")

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

        # Create per-sample working directory
        mkdir -p {params.workdir}
        TMPDIR="{params.workdir}/tmp"
        mkdir -p "$TMPDIR"
        SEQDB="{params.workdir}/seqDB"
        CLUDB="{params.workdir}/cluDB"
        REPDB="{params.workdir}/repDB"
        
        echo "[$(date)] START dedup: {wildcards.sample}" > {log}
        
 
                
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
        OUTPUT_COUNT=$(zcat {output.fasta} | grep -c ">" {output.fasta} || true)
        # REMOVED=$((SEQ_COUNT - OUTPUT_COUNT))
        # echo "Input sequences  : $SEQ_COUNT"  2>> {log}
        echo "Unique sequences : $OUTPUT_COUNT" 2>> {log}
        # echo "Removed          : $REMOVED"      2>> {log}
        
        """


