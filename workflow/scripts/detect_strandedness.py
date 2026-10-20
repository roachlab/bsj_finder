#!/usr/bin/env python3
"""
Per-sample library strandedness, from the reads and the reference alone.

METHOD. A read spanning a splice junction pins the transcript's orientation via the
GT-AG motif, and minimap2 (-x splice:sr) reports that as ts:A, RELATIVE TO THE READ.
So the ts:A distribution over R1 alone answers the question:

    R1 mostly ts:-   -> R1 is antisense  -> reverse  (dUTP / TruSeq Stranded)
    R1 mostly ts:+   -> R1 is sense      -> forward  (ligation-based)
    roughly 50/50    -> unstranded

No annotation is needed -- only the genome already being aligned against. R1 alone is
used deliberately: R2 is reverse-complemented in rename_read_ids_pe, which erases the
mate signature, and using R1 makes the test identical for single-end samples.

VALIDATION. On GBM sample 150 this gives R1 96.2% antisense. Scoring exonic reads
against GENCODE v44 gene strand independently gives 97.2% antisense -- agreement
within a percentage point. That cross-check is what establishes ts:A is read-relative
(a reference-relative tag would sit near 50/50, since transcripts occur on both
genomic strands).

SAMPLE SIZE. The verdict is stable from ~5k reads: GBM R1 reads 5k/20k/50k/150k/300k
give 97.1/97.0/96.3/96.2/96.2% ts:-. The default of 50k is well past saturation and
costs a second or two against a prebuilt index.

IT REFUSES TO GUESS. Too few spliced reads, or a split that is neither dominant nor
balanced, is a hard error rather than a default. A wrong call is invisible and puts
every junction on the wrong strand, so failing the run is the safer outcome -- an
RNase-R-treated circRNA library with little linear splicing is exactly the case that
should stop rather than be guessed at.
"""
import gzip
import subprocess
import sys


def tabulate(mmi, reads, n_reads, threads, log):
    """Align a subsample of R1 and count ts:A over primary alignments."""
    fq = []
    with gzip.open(reads, "rt") as fh:
        for i, line in enumerate(fh):
            if i >= n_reads * 4:
                break
            fq.append(line)
    if not fq:
        return 0, 0, 0
    n_sampled = len(fq) // 4

    proc = subprocess.Popen(
        ["minimap2", "-ax", "splice:sr", "-t", str(threads), "--secondary=no", mmi, "-"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log, text=True,
    )
    plus = minus = 0
    try:
        proc.stdin.write("".join(fq))
        proc.stdin.close()
        for line in proc.stdout:
            if line[0] == "@":
                continue
            f = line.split("\t")
            flag = int(f[1])
            if flag & 0x904:                      # unmapped / secondary / supplementary
                continue
            for tag in f[11:]:
                if tag.startswith("ts:A:"):
                    if tag[5] == "+":
                        plus += 1
                    elif tag[5] == "-":
                        minus += 1
                    break
    finally:
        proc.stdout.close()
        rc = proc.wait()
    if rc != 0:
        raise SystemExit(f"minimap2 failed (exit {rc}) during strandedness detection")
    return n_sampled, plus, minus


def call(plus, minus, min_spliced, min_frac):
    """(verdict, reason). Ambiguity is an error, never a default."""
    total = plus + minus
    if total < min_spliced:
        return None, (f"only {total} spliced reads with a ts:A tag (need {min_spliced}). "
                      "Too few linear splice junctions to call strandedness -- raise "
                      "strandedness_reads, or declare library_strandedness explicitly.")
    f_minus = minus / total
    if f_minus >= min_frac:
        return "reverse", f"R1 {100*f_minus:.1f}% antisense"
    if (1 - f_minus) >= min_frac:
        return "forward", f"R1 {100*(1-f_minus):.1f}% sense"
    if 0.35 <= f_minus <= 0.65:
        return "unstranded", f"R1 {100*f_minus:.1f}% antisense -- balanced"
    return None, (f"ambiguous: R1 {100*f_minus:.1f}% antisense, which is neither "
                  f"dominant (>={100*min_frac:.0f}%) nor balanced (35-65%). Declare "
                  "library_strandedness explicitly if you know the protocol.")


def main(sample, reads, mmi, out, declared, n_reads, min_spliced, min_frac, threads, logpath):
    with open(logpath, "w") as log:
        n_sampled, plus, minus = tabulate(mmi, reads, n_reads, threads, log)
        verdict, reason = call(plus, minus, min_spliced, min_frac)
        total = plus + minus
        log.write(f"[strandedness] {sample}: sampled {n_sampled} R1 reads, "
                  f"{total} with ts:A (+{plus} / -{minus})\n")
        log.write(f"[strandedness] verdict: {verdict or 'UNDETERMINED'} -- {reason}\n")
        if verdict is None:
            raise SystemExit(f"{sample}: could not determine strandedness -- {reason}")
        if declared != "auto" and declared != verdict:
            # The config asserted something the reads contradict. That is worth stopping
            # for: one of the two is wrong, and guessing which would be worse than either.
            msg = (f"{sample}: library_strandedness is declared '{declared}' but the reads "
                   f"say '{verdict}' ({reason}). Fix the config, or set it to 'auto'.")
            log.write(f"[strandedness] MISMATCH: {msg}\n")
            raise SystemExit(msg)
        with open(out, "w") as fo:
            fo.write("sample\tverdict\tdeclared\treads_sampled\tspliced_ts\tts_plus\tts_minus\n")
            fo.write(f"{sample}\t{verdict}\t{declared}\t{n_sampled}\t{total}\t{plus}\t{minus}\n")


if "snakemake" in globals():
    smk = globals()["snakemake"]                      # noqa: F821  (injected by script:)
    main(
        sample=smk.wildcards.sample,
        reads=str(smk.input.r1),
        mmi=str(smk.input.mmi),
        out=str(smk.output.tsv),
        declared=str(smk.params.declared),
        n_reads=int(smk.params.n_reads),
        min_spliced=int(smk.params.min_spliced),
        min_frac=float(smk.params.min_frac),
        threads=int(smk.threads),
        logpath=str(smk.log[0]),
    )
elif __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Detect library strandedness from reads + reference")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--reads", required=True, help="R1 fastq.gz (R1 is never reverse-complemented)")
    ap.add_argument("--mmi", required=True, help="minimap2 index (or a fasta)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--log", default="/dev/stderr")
    ap.add_argument("--declared", default="auto", choices=("auto", "forward", "reverse", "unstranded"))
    ap.add_argument("--n-reads", type=int, default=50000)
    ap.add_argument("--min-spliced", type=int, default=200)
    ap.add_argument("--min-frac", type=float, default=0.80)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    main(a.sample, a.reads, a.mmi, a.out, a.declared, a.n_reads,
         a.min_spliced, a.min_frac, a.threads, a.log)
