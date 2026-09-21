# Read trimming (fastp) and read-ID renaming (seqtk).
#
# Both steps are split into _pe/_se variants. A single rule cannot declare a
# different NUMBER of outputs per sample, so serving both layouts from one rule
# meant R2/RS had to be smuggled through as params -- untracked by Snakemake,
# therefore unreachable by temp() and never cleaned up. The split makes every
# file a real output. PE_SAMPLES/SE_SAMPLES and sample_regex() come from the
# Snakefile; wildcard_constraints keep the two variants from overlapping.
#
# Every output here is temp(): both the fastp files and the renamed copies are
# intermediates, and on a full-size run they dwarf the results that are actually
# kept (~200 GB for the 38-sample PRJEB34650 set). The renamed reads are freed once
# polished_bsj_seq -- their last consumer -- has run for that sample; the fastp
# files go as soon as the rename rule has. Because they are temp, re-running any
# downstream step later re-trims from the source FASTQs. Pass --notemp to keep
# them for a run where you want to inspect the trimmed reads.
#
# fastp flags are identical across the two beyond the pairing-specific ones
# (-I/-O/--unpaired*/--detect_adapter_for_pe). Keep them in sync.
#
# --length_required: see config key fastp_min_len. It MUST sit below the library's
# read length -- above it, fastp discards every read, still exits 0, and leaves
# valid-but-empty FASTQs that the rest of the DAG happily consumes (this silently
# emptied all 37 paired 75 bp PRJEB34650 samples under the old hardcoded 90).
# --cut_tail shortens reads further, so leave headroom below the raw read length.


rule trimnami_fastp_pe:
    input:
        r1 = lambda wc: samples[wc.sample]["R1"],
        r2 = lambda wc: samples[wc.sample]["R2"],
    output:
        r1     = temp(os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz")),
        r2     = temp(os.path.join(out_dir, "trim_results", "{sample}_R2_trim.fastq.gz")),
        orphan = temp(os.path.join(out_dir, "trim_results", "{sample}_RS_trim.fastq.gz")),
    wildcard_constraints:
        sample = sample_regex(PE_SAMPLES),
    params:
        min_len = config.get("fastp_min_len", 90),
    conda:
        os.path.join("..", "env", "trimnami.yaml")
    log:
        os.path.join("logs", "trimnami_fastp.{sample}.log")
    threads:
        8
    shell:
        """
        echo "[$(date)] Running PAIRED-END trimming for {wildcards.sample}" >> {log}

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --unpaired1 {output.orphan} \
            --unpaired2 {output.orphan} \
            -z 1 \
            --thread {threads} \
            --qualified_quality_phred 15 \
            --length_required {params.min_len} \
            --cut_tail \
            --cut_tail_window_size 25 \
            --cut_tail_mean_quality 15 \
            --dedup \
            --dup_calc_accuracy 4 \
            --trim_poly_x \
            --detect_adapter_for_pe \
            2>> {log}

        # fastp only creates the --unpaired file if at least one read is orphaned.
        # It is a declared output now, so guarantee it exists.
        #
        # NB: this placeholder is a ~20-byte valid gzip of nothing, so it is NOT
        # empty -- the downstream `[ -s ]` guards (minimap2_alig, polished_bsj_seq)
        # see it as present and pass it through. That is harmless, because
        # decompressing an empty member yields no reads, but `-s` here means "the
        # file exists", NOT "the file has reads". Do not write a new guard on the
        # assumption that it means the latter.
        [ -f "{output.orphan}" ] || gzip -1 < /dev/null > "{output.orphan}"
        """


rule trimnami_fastp_se:
    input:
        r1 = lambda wc: samples[wc.sample]["R1"],
    output:
        r1 = temp(os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz")),
    wildcard_constraints:
        sample = sample_regex(SE_SAMPLES),
    params:
        min_len = config.get("fastp_min_len", 90),
    conda:
        os.path.join("..", "env", "trimnami.yaml")
    log:
        os.path.join("logs", "trimnami_fastp.{sample}.log")
    threads:
        8
    shell:
        """
        echo "[$(date)] Running SINGLE-END trimming for {wildcards.sample}" >> {log}

        fastp \
            -i {input.r1} \
            -o {output.r1} \
            -z 1 \
            --thread {threads} \
            --qualified_quality_phred 15 \
            --length_required {params.min_len} \
            --cut_tail \
            --cut_tail_window_size 25 \
            --cut_tail_mean_quality 15 \
            --dedup \
            --dup_calc_accuracy 4 \
            --trim_poly_x \
            2>> {log}
        """


# seqtk rename rewrites read IDs to {sample}.{mate}.{n} so reads stay attributable
# once all three files are concatenated into one single-end stream for alignment.

rule rename_read_ids_pe:
    input:
        r1     = os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz"),
        r2     = os.path.join(out_dir, "trim_results", "{sample}_R2_trim.fastq.gz"),
        orphan = os.path.join(out_dir, "trim_results", "{sample}_RS_trim.fastq.gz"),
    output:
        r1     = temp(os.path.join(out_dir, "trim_results", "{sample}_R1_trim_rename.fastq.gz")),
        r2     = temp(os.path.join(out_dir, "trim_results", "{sample}_R2_trim_rename.fastq.gz")),
        orphan = temp(os.path.join(out_dir, "trim_results", "{sample}_RS_trim_rename.fastq.gz")),
    wildcard_constraints:
        sample = sample_regex(PE_SAMPLES),
    conda:
        os.path.join("..", "env", "seqtk.yaml")
    log:
        os.path.join("logs", "rename_read_ids.{sample}.log")
    shell:
        """
        seqtk rename {input.r1}     '{wildcards.sample}.1.' | gzip -1 - > {output.r1}
        seqtk rename {input.r2}     '{wildcards.sample}.2.' | gzip -1 - > {output.r2}
        seqtk rename {input.orphan} '{wildcards.sample}.S.' | gzip -1 - > {output.orphan}
        """


rule rename_read_ids_se:
    input:
        r1 = os.path.join(out_dir, "trim_results", "{sample}_R1_trim.fastq.gz"),
    output:
        r1 = temp(os.path.join(out_dir, "trim_results", "{sample}_R1_trim_rename.fastq.gz")),
    wildcard_constraints:
        sample = sample_regex(SE_SAMPLES),
    conda:
        os.path.join("..", "env", "seqtk.yaml")
    log:
        os.path.join("logs", "rename_read_ids.{sample}.log")
    shell:
        """
        seqtk rename {input.r1} '{wildcards.sample}.1.' | gzip -1 - > {output.r1}
        """
