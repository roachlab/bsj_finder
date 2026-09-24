# ═══════════════════════════════════════════════════════════════════════════
# Pool dedup — global BSJ pool across ALL samples
#
# concat_all_bsj                  : concatenate the per-sample raw BSJ fastas
#                                   (gzip) into one multi-member gzip
# mmseqs2_dedup_pool              : mmseqs2 linclust over the combined RAW pool
# concat_all_polished             : concatenate the per-sample POLISHED fastas
#                                   (plain text) into one file
# mmseqs2_dedup_pool_polished     : mmseqs2 linclust over the combined POLISHED
#                                   pool — the final cross-sample unique set of
#                                   polished BSJs
#
# Same mmseqs params as the per-sample rule (mmseqs_min_seq_id / min_cov /
# cov_mode) and the same conda-env-or-static-binary selection (see the MMSEQS2
# BINARY SELECTION block in the Snakefile).
# ═══════════════════════════════════════════════════════════════════════════

rule concat_all_bsj:
    input:
        lambda wildcards: expand(
            os.path.join(out_dir, "bsj_results", "{sample}_bsj_sequences.fasta.gz"),
            sample=SAMPLES
        )
    output:
        os.path.join(out_dir, "pool_results", "all_samples.bsj_sequences.fasta.gz")
    log:
        os.path.join("logs", "concat_all_bsj.log")
    shell:
        """
        # Per-sample fastas are gzip; concatenating gzip members yields a
        # valid multi-member gzip file that zcat/gunzip handle transparently.
        cat {input} > {output}
        echo "[$(date)] concatenated {input} -> {output}" > {log}
        echo "[$(date)] total sequences: $(zcat {output} | grep -c '^>')" >> {log}
        """


rule mmseqs2_dedup_pool:
    input:
        fasta = os.path.join(out_dir, "pool_results", "all_samples.bsj_sequences.fasta.gz")
    output:
        fasta = os.path.join(out_dir, "pool_results", "all_samples.bsj_pool_dedup.fasta"),
        tsv   = os.path.join(out_dir, "pool_results", "all_samples.bsj_pool_cluster.tsv")

    # mmseqs comes from the conda env by default; with `mmseqs_static: true` the
    # env is dropped (conda: None) and the vendored SSE4.1 binary is used instead --
    # see the MMSEQS2 BINARY SELECTION block in the Snakefile.
    conda: MMSEQS_CONDA

    # threads: 1 is a CORRECTNESS setting, not a performance one.
    #
    # mmseqs linclust is order-dependent when multithreaded: two identical runs pick
    # different cluster representatives. Measured on this data, rep counts moved by
    # ~0.03-0.07% run to run (GBM polished pool 763,313 vs 763,542), which propagates
    # into polished_results and every n_samples in the summary table -- so no published
    # junction catalogue could be reproduced exactly. Single-threaded linclust is
    # deterministic; verified by two full end-to-end runs producing a byte-identical
    # summary table.
    #
    # The cost is ~12 min per dataset, all of it in the two pool rules (measured on the
    # GBM pools: polished 5m09s vs 0m43s at 16 threads, raw 9m07s vs 1m26s). The
    # per-sample rule costs almost nothing, since --cores still runs many samples at
    # once. Peak scratch is slightly LOWER single-threaded, so there is no memory
    # argument against it either.
    #
    # If you ever trade this back for speed, expect the catalogue to shift slightly and
    # say so wherever the numbers are published.
    threads:
        1

    log:
        os.path.join("logs", "mmseqs2_dedup_pool.log")

    params:
        min_seq_id = config.get("mmseqs_min_seq_id", 0.95),
        min_cov = config.get("mmseqs_min_cov",    0.70),
        cov_mode = config.get("mmseqs_cov_mode",   0),
        workdir = os.path.join(out_dir, "pool_results", "all_samples_mmseqs2_tmp"),
        mmseqs_bin_dir = MMSEQS_BIN_DIR,

    shell:
        """
        set -euo pipefail

        # The static SSE4.1 build, when enabled, takes precedence over any mmseqs
        # already on PATH (including the conda env's). Empty => leave PATH alone.
        if [[ -n "{params.mmseqs_bin_dir}" ]]; then
            export PATH="{params.mmseqs_bin_dir}:$PATH"
        fi

        mkdir -p {params.workdir}
        TMPDIR="{params.workdir}/tmp"
        mkdir -p "$TMPDIR"
        SEQDB="{params.workdir}/seqDB"
        CLUDB="{params.workdir}/cluDB"
        REPDB="{params.workdir}/repDB"

        echo "[$(date)] START pool dedup (all samples, raw)" > {log}

        SEQ_COUNT=$(gunzip -c {input.fasta} | grep -c "^>" || echo 0)
        SEQ_COUNT=$(echo "$SEQ_COUNT" | tr -d '[:space:]')
        echo "[$(date)] Input sequence count: $SEQ_COUNT" >> {log}

        if [[ "$SEQ_COUNT" -eq 0 ]]; then
            echo "[$(date)] No sequences found — writing empty outputs" >> {log}
            gunzip -c {input.fasta} > {output.fasta}
            printf "query\\ttarget\\n" > {output.tsv}
            rm -rf {params.workdir}
            exit 0
        fi

        if [[ "$SEQ_COUNT" -eq 1 ]]; then
            echo "[$(date)] Only 1 sequence — skipping clustering" >> {log}
            gunzip -c {input.fasta} > {output.fasta}
            gunzip -c {input.fasta} | grep "^>" | sed 's/^>//' | \\
                awk '{{print $1"\\t"$1}}' > {output.tsv}
            rm -rf {params.workdir}
            exit 0
        fi

        mmseqs createdb \\
            {input.fasta} \\
            "$SEQDB" \\
            --threads {threads} \\
            2>> {log}

        mmseqs linclust \\
            "$SEQDB" \\
            "$CLUDB" \\
            "$TMPDIR" \\
            --min-seq-id {params.min_seq_id} \\
            -c {params.min_cov}    \\
            --cov-mode {params.cov_mode}   \\
            --threads {threads} \\
            2>> {log}

        mmseqs result2repseq \\
            "$SEQDB" \\
            "$CLUDB" \\
            "$REPDB" \\
            --threads {threads} \\
            2>> {log}

        mmseqs convert2fasta \\
            "$REPDB" \\
            {output.fasta} \\
            2>> {log}

        mmseqs createtsv \\
            "$SEQDB" \\
            "$SEQDB" \\
            "$CLUDB" \\
            {output.tsv} \\
            --threads {threads} \\
            2>> {log}

        rm -rf {params.workdir}

        # Summary
        OUTPUT_COUNT=$(grep -c ">" {output.fasta} || echo 0)
        echo "[$(date)] Input sequences  : $SEQ_COUNT" >> {log}
        echo "[$(date)] Unique sequences : $OUTPUT_COUNT" >> {log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# Summary table — single BSJ table across all samples
#
# Inner join: the final unique set of polished BSJs (representatives from the
# cross-sample pool dedup) is joined back onto the per-sample bsj_cor rows, so
# only junctions that SURVIVED deduplication appear. Sample of origin is dropped
# — the same junction called in several samples collapses to one row; n_samples
# counts the distinct samples contributing members to its cluster(s).
# Stdlib-only script (no conda env).
# ─────────────────────────────────────────────────────────────────────────────

rule bsj_summary_table:
    input:
        cor = lambda wildcards: expand(
            os.path.join(out_dir, "bsj_results", "{sample}.bsj_cor.tsv.gz"),
            sample=SAMPLES
        ),
        pool_fasta = os.path.join(out_dir, "pool_results", "all_samples.polished_pool_dedup.fasta"),
        pool_cluster = os.path.join(out_dir, "pool_results", "all_samples.polished_pool_cluster.tsv")
    output:
        tsv = os.path.join(out_dir, "pool_results", "all_samples.bsj_summary.tsv")
    log:
        os.path.join("logs", "bsj_summary_table.log")
    params:
        locus_tol = config.get("bsj_locus_tol", 10),
    # Stdlib-only script (like sa_detect.py) -- no conda env required, so it runs under
    # Snakemake's own interpreter.
    #
    # script: rather than shell: -- the script then reads input/output/params/log off
    # the injected `snakemake` object instead of a hand-built argv, and, importantly,
    # Snakemake HASHES THE FILE as part of the rule's code. Under the old
    # `python3 {params.script}` form the script was just a param string, so edits to it
    # never tripped the code rerun trigger and a changed summary script would silently
    # not re-run on an otherwise up-to-date tree.
    #
    # The path is resolved relative to THIS .smk file (same convention as conda:).
    script:
        os.path.join("..", "scripts", "bsj_summary.py")


# ─────────────────────────────────────────────────────────────────────────────
# POLISHED pool — cross-sample unique polished BSJs (final step)
# ─────────────────────────────────────────────────────────────────────────────

rule concat_all_polished:
    input:
        lambda wildcards: expand(
            os.path.join(out_dir, "polished_results", "{sample}.polished.fasta"),
            sample=SAMPLES
        )
    output:
        os.path.join(out_dir, "pool_results", "all_samples.polished.fasta")
    log:
        os.path.join("logs", "concat_all_polished.log")
    shell:
        """
        # Per-sample polished fastas are plain text; a straight cat is a valid
        # combined fasta. Samples with zero polished sequences are empty files
        # and contribute nothing.
        cat {input} > {output}
        echo "[$(date)] concatenated {input} -> {output}" > {log}
        echo "[$(date)] total sequences: $(grep -c '^>' {output})" >> {log}
        """


rule mmseqs2_dedup_pool_polished:
    input:
        fasta = os.path.join(out_dir, "pool_results", "all_samples.polished.fasta")
    output:
        fasta = os.path.join(out_dir, "pool_results", "all_samples.polished_pool_dedup.fasta"),
        tsv   = os.path.join(out_dir, "pool_results", "all_samples.polished_pool_cluster.tsv")

    # mmseqs comes from the conda env by default; with `mmseqs_static: true` the
    # env is dropped (conda: None) and the vendored SSE4.1 binary is used instead --
    # see the MMSEQS2 BINARY SELECTION block in the Snakefile.
    conda: MMSEQS_CONDA

    # threads: 1 is a CORRECTNESS setting, not a performance one.
    #
    # mmseqs linclust is order-dependent when multithreaded: two identical runs pick
    # different cluster representatives. Measured on this data, rep counts moved by
    # ~0.03-0.07% run to run (GBM polished pool 763,313 vs 763,542), which propagates
    # into polished_results and every n_samples in the summary table -- so no published
    # junction catalogue could be reproduced exactly. Single-threaded linclust is
    # deterministic; verified by two full end-to-end runs producing a byte-identical
    # summary table.
    #
    # The cost is ~12 min per dataset, all of it in the two pool rules (measured on the
    # GBM pools: polished 5m09s vs 0m43s at 16 threads, raw 9m07s vs 1m26s). The
    # per-sample rule costs almost nothing, since --cores still runs many samples at
    # once. Peak scratch is slightly LOWER single-threaded, so there is no memory
    # argument against it either.
    #
    # If you ever trade this back for speed, expect the catalogue to shift slightly and
    # say so wherever the numbers are published.
    threads:
        1

    log:
        os.path.join("logs", "mmseqs2_dedup_pool_polished.log")

    params:
        min_seq_id = config.get("mmseqs_min_seq_id", 0.95),
        min_cov = config.get("mmseqs_min_cov",    0.70),
        cov_mode = config.get("mmseqs_cov_mode",   0),
        workdir = os.path.join(out_dir, "pool_results", "all_samples_polished_mmseqs2_tmp"),
        mmseqs_bin_dir = MMSEQS_BIN_DIR,

    shell:
        """
        set -euo pipefail

        # The static SSE4.1 build, when enabled, takes precedence over any mmseqs
        # already on PATH (including the conda env's). Empty => leave PATH alone.
        if [[ -n "{params.mmseqs_bin_dir}" ]]; then
            export PATH="{params.mmseqs_bin_dir}:$PATH"
        fi

        mkdir -p {params.workdir}
        TMPDIR="{params.workdir}/tmp"
        mkdir -p "$TMPDIR"
        SEQDB="{params.workdir}/seqDB"
        CLUDB="{params.workdir}/cluDB"
        REPDB="{params.workdir}/repDB"

        echo "[$(date)] START polished pool dedup (all samples)" > {log}

        SEQ_COUNT=$(grep -c "^>" {input.fasta} || echo 0)
        SEQ_COUNT=$(echo "$SEQ_COUNT" | tr -d '[:space:]')
        echo "[$(date)] Input sequence count: $SEQ_COUNT" >> {log}

        if [[ "$SEQ_COUNT" -eq 0 ]]; then
            echo "[$(date)] No sequences found — writing empty outputs" >> {log}
            touch {output.fasta}
            printf "query\\ttarget\\n" > {output.tsv}
            rm -rf {params.workdir}
            exit 0
        fi

        if [[ "$SEQ_COUNT" -eq 1 ]]; then
            echo "[$(date)] Only 1 sequence — skipping clustering" >> {log}
            cp {input.fasta} {output.fasta}
            grep "^>" {input.fasta} | sed 's/^>//' | \\
                awk '{{print $1"\\t"$1}}' > {output.tsv}
            rm -rf {params.workdir}
            exit 0
        fi

        mmseqs createdb \\
            {input.fasta} \\
            "$SEQDB" \\
            --threads {threads} \\
            2>> {log}

        mmseqs linclust \\
            "$SEQDB" \\
            "$CLUDB" \\
            "$TMPDIR" \\
            --min-seq-id {params.min_seq_id} \\
            -c {params.min_cov}    \\
            --cov-mode {params.cov_mode}   \\
            --threads {threads} \\
            2>> {log}

        mmseqs result2repseq \\
            "$SEQDB" \\
            "$CLUDB" \\
            "$REPDB" \\
            --threads {threads} \\
            2>> {log}

        mmseqs convert2fasta \\
            "$REPDB" \\
            {output.fasta} \\
            2>> {log}

        mmseqs createtsv \\
            "$SEQDB" \\
            "$SEQDB" \\
            "$CLUDB" \\
            {output.tsv} \\
            --threads {threads} \\
            2>> {log}

        rm -rf {params.workdir}

        # Summary
        OUTPUT_COUNT=$(grep -c ">" {output.fasta} || echo 0)
        echo "[$(date)] Input sequences  : $SEQ_COUNT" >> {log}
        echo "[$(date)] Unique sequences : $OUTPUT_COUNT" >> {log}
        """
