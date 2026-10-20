# Shared minimap2 index + per-sample strandedness detection.
#
# The index rule exists for both reasons: it makes strandedness detection nearly free
# (otherwise every sample would rebuild a multi-Gb index for a 50k-read alignment), and
# it removes the per-sample index build from minimap2_alig, which previously rebuilt it
# 38 times on the Tasmanian devil set.
#
# The index is built with the config's -k/-w. minimap2 then takes those parameters FROM
# the index and warns that -k/-w are overridden -- which is correct and expected, since
# both consumers want the same seeding.

rule minimap2_index:
    input:
        genome = config['genome']
    output:
        mmi = os.path.join(out_dir, "ref", "genome.mmi")
    conda:
        os.path.join("..", "env", "minimap2.yaml")
    log:
        os.path.join("logs", "minimap2_index.log")
    threads:
        8
    params:
        k = config.get("minimap2_k", 11),
        w = config.get("minimap2_w", 5),
    shell:
        """
        mkdir -p "$(dirname {output.mmi})"
        minimap2 -x sr -k {params.k} -w {params.w} -t {threads} \
            -d {output.mmi} {input.genome} 2> {log}
        """


# Per-sample, not per-dataset: a single accession list in this project turned out to
# contain two species AND two protocols, so a dataset-level assumption is not safe.
#
# Uses R1 only. R2 is reverse-complemented in rename_read_ids_pe, which erases the mate
# signature, and R1-only makes the test identical for single-end samples.
rule detect_strandedness:
    input:
        r1  = os.path.join(out_dir, "trim_results", "{sample}_R1_trim_rename.fastq.gz"),
        mmi = os.path.join(out_dir, "ref", "genome.mmi"),
    output:
        tsv = os.path.join(out_dir, "qc_results", "{sample}.strandedness.tsv")
    conda:
        os.path.join("..", "env", "minimap2.yaml")
    log:
        os.path.join("logs", "detect_strandedness.{sample}.log")
    threads:
        4
    params:
        declared    = STRANDEDNESS,                                  # validated in the Snakefile
        n_reads     = config.get("strandedness_reads",       50000), # verdict saturates by ~5k
        min_spliced = config.get("strandedness_min_spliced",   200), # below this: refuse to call
        min_frac    = config.get("strandedness_min_frac",      0.80),
    script:
        os.path.join("..", "scripts", "detect_strandedness.py")
