#!/usr/bin/env python3
"""
Chimeric back-splice-junction detector for minimap2 SAM output (SA-tag / soft-clip based).

Motivation: the PAF + `-x splice:sr` path surfaces only ~half of true junction reads as
two-arm splits. Running minimap2 in SAM mode (`-a`) and reading its `SA:Z` supplementary
alignments (the same signal CIRI2 reads from BWA-MEM) recovers ~3x more arms while keeping
minimap2 as the aligner (speed + long-read capability).

For each read it pairs the primary alignment with each supplementary segment on the same
chromosome and strand, checks the "chiastic" out-of-order signature of a back-splice, and
emits one row per supporting read in bsj_finder's existing bsj_cor.tsv schema so the current
converter (bsj_finder_to_bed.py) and benchmark score it unchanged:

    q_id  bsj_coord  chromosome  s_start  s_end  s_strand  overlap_trimmed_nt

Coordinates are 0-based, s_start<s_end (circle_start / circle_end convention of truth_bsj.tsv).
Dependency-free (stdlib only). Reads SAM from a file arg or stdin.
"""
import sys, re, argparse

_CIG = re.compile(r'(\d+)([MIDNSHP=X])')

def ref_len(cig):
    return sum(int(n) for n, op in _CIG.findall(cig) if op in "MDN=X")

def query_aln_len(cig):
    return sum(int(n) for n, op in _CIG.findall(cig) if op in "MI=X")

def clips(cig):
    """leading, trailing soft/hard clip lengths (as written on the aligned strand)."""
    lead = 0
    m = re.match(r'^(\d+)([SH])', cig)
    if m: lead = int(m.group(1))
    trail = 0
    m = re.search(r'(\d+)([SH])$', cig)
    if m: trail = int(m.group(1))
    return lead, trail

def seg_from(chrom, pos1, strand, cig, mapq):
    """Build a segment dict in ORIGINAL-read coordinates.
    ref span is 0-based half-open [ref_start, ref_end)."""
    rs = pos1 - 1
    re_ = rs + ref_len(cig)
    lead, trail = clips(cig)
    qlen = query_aln_len(cig)
    # On '-' strand the CIGAR is written for the revcomp read, so the leading clip
    # corresponds to the 3' end of the original read: flip to original-read frame.
    if strand == '+':
        q_start = lead
    else:
        q_start = trail
    q_end = q_start + qlen
    return dict(chrom=chrom, strand=strand, rs=rs, re=re_, qs=q_start, qe=q_end, mapq=mapq)

def detect(sam_fh, min_anchor=30, min_mapq=1, max_gap_size=100000,
           max_read_gap=15, max_read_ovlp=15):
    for ln in sam_fh:
        if ln.startswith("@"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0x4:                      # unmapped
            continue
        if flag & 0x100 or flag & 0x800:    # only anchor on the primary line
            continue
        qname, rname, pos, mapq, cig = f[0], f[2], int(f[3]), int(f[4]), f[5]
        pstrand = '-' if (flag & 0x10) else '+'
        sa = None
        for t in f[11:]:
            if t.startswith("SA:Z:"):
                sa = t[5:]; break
        if not sa:
            continue
        prim = seg_from(rname, pos, pstrand, cig, mapq)
        for entry in sa.split(";"):
            if not entry:
                continue
            p = entry.split(",")
            if len(p) < 6:
                continue
            srname, spos, sstrand, scig, smapq, _nm = p[0], int(p[1]), p[2], p[3], int(p[4]), p[5]
            if srname != prim["chrom"] or sstrand != prim["strand"]:
                continue                    # different chr / inversion -> not a standard BSJ
            supp = seg_from(srname, spos, sstrand, scig, smapq)

            # order the two arms by their position in the ORIGINAL read
            a, b = (prim, supp) if prim["qs"] <= supp["qs"] else (supp, prim)
            # anchor length & mapq filters on both arms
            if (a["qe"] - a["qs"]) < min_anchor or (b["qe"] - b["qs"]) < min_anchor:
                continue
            if a["mapq"] < min_mapq or b["mapq"] < min_mapq:
                continue
            # read-space gap/overlap between arms at the junction
            read_gap = a["qe"] - b["qs"]     # >0 overlap (microhomology), <0 gap
            if read_gap > max_read_ovlp or read_gap < -max_read_gap:
                continue

            strand = prim["strand"]
            # inner (junction-facing) genomic coordinate of each arm
            if strand == '+':
                # arm 'a' (first in read) = donor side, inner = its ref_end
                # arm 'b' (second in read) = acceptor side, inner = its ref_start
                bp_donor = a["re"]           # 0-based half-open end
                bp_accept = b["rs"]          # 0-based start
                # back-splice requires out-of-order: donor downstream of acceptor
                if not (a["rs"] > b["rs"]):
                    continue
                left, right = bp_accept, bp_donor
            else:
                # minus strand: geometry mirrors
                bp_donor = a["rs"]
                bp_accept = b["re"]
                if not (a["rs"] < b["rs"]):
                    continue
                left, right = bp_donor, bp_accept
            if right <= left:
                continue
            if (right - left) > max_gap_size:
                continue
            # Microhomology: when the two arms overlap by `ovl` bases in the read, those
            # bases are aligned by both arms and the exact breakpoint is ambiguous over
            # that window. minimap assigns the overlap to the downstream (right) arm, so
            # trim it from `right` to match the dataset's left-anchored breakpoint
            # convention (keeps canonical calls exact, not +ovl-shifted).
            ovl = read_gap if read_gap > 0 else 0
            if ovl:
                right = right - ovl
            if right <= left:
                continue
            # read-space position of the junction (end of the first arm), used by the
            # downstream extract_bsj_seq step to window the read around the BSJ.
            bsj_read_pos = max(0, a["qe"])
            yield f"{qname}\t{bsj_read_pos}\t{prim['chrom']}\t{left}\t{right}\t{strand}\t{ovl}\n"

HEADER = "q_id\tbsj_coord\tchromosome\ts_start\ts_end\ts_strand\toverlap_trimmed_nt\n"


def write_table(sam_fh, out_fh, min_anchor, min_mapq, max_gap_size,
                max_read_gap, max_read_ovlp):
    """Stream SAM records in, bsj_cor rows out. Returns the row count."""
    out_fh.write(HEADER)
    n = 0
    for line in detect(sam_fh, min_anchor, min_mapq, max_gap_size,
                       max_read_gap, max_read_ovlp):
        out_fh.write(line)
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="SA-tag BSJ detector for minimap2 SAM")
    ap.add_argument("--sam", default="-", help="SAM file, or '-' for stdin")
    ap.add_argument("--out", default="-", help="output bsj_cor.tsv, or '-' for stdout")
    ap.add_argument("--min-anchor", type=int, default=30)
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--max-gap-size", type=int, default=100000)
    ap.add_argument("--max-read-gap", type=int, default=15)
    ap.add_argument("--max-read-ovlp", type=int, default=15)
    a = ap.parse_args()
    fh = sys.stdin if a.sam == "-" else open(a.sam)
    out = sys.stdout if a.out == "-" else open(a.out, "w")
    n = 0
    try:
        n = write_table(fh, out, a.min_anchor, a.min_mapq, a.max_gap_size,
                        a.max_read_gap, a.max_read_ovlp)
    finally:
        if fh is not sys.stdin: fh.close()
        if out is not sys.stdout: out.close()
    sys.stderr.write(f"[sa_detect] {n} BSJ read-support rows\n")


def _from_snakemake(smk):
    """Run under the `script:` directive.

    The rule used to be a shell pipeline (samtools | python | gzip). The samtools and
    gzip ends move in here so the file itself is the rule's code -- Snakemake then
    hashes it, and edits to the detector actually trip the rerun trigger. Streaming is
    preserved: samtools stays a subprocess whose stdout is consumed lazily, so the SAM
    is never held in memory.

    -F 0x900 drops secondary/supplementary records. detect() anchors on the primary
    line and rebuilds the second arm from its SA:Z string, so those records are ~half
    the BAM and 100% wasted parsing -- but they ARE the second arm, so they stay in the
    BAM on disk for IGV. samtools comes from the rule's conda env.
    """
    import gzip, os, subprocess

    p = smk.params
    out_path = str(smk.output.tsv)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    with open(str(smk.log[0]), "w") as log_fh:
        cmd = ["samtools", "view", "-@", str(smk.threads), "-F", "0x900",
               str(smk.input.bam)]
        log_fh.write(f"[sa_detect] {' '.join(cmd)}\n")
        log_fh.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log_fh, text=True)
        try:
            with gzip.open(out_path, "wt") as out_fh:
                n = write_table(proc.stdout, out_fh,
                                p.min_anchor, p.min_mapq, p.max_gap_size,
                                p.max_read_gap, p.max_read_ovlp)
        finally:
            proc.stdout.close()
            rc = proc.wait()
        if rc != 0:
            # Without this the rule would "succeed" on a truncated BAM: the pipeline's
            # old `samtools | python` form hid samtools' exit status the same way.
            raise SystemExit(f"samtools view failed (exit {rc}) -- see {smk.log[0]}")
        log_fh.write(f"[sa_detect] {n} BSJ read-support rows\n")


if "snakemake" in globals():
    _from_snakemake(globals()["snakemake"])   # noqa: F821  (injected by script:)
elif __name__ == "__main__":
    main()
