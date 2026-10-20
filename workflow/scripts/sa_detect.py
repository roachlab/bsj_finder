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

    q_id  bsj_coord  chromosome  s_start  s_end  s_strand  overlap_trimmed_nt  mapq

`mapq` is the MAPQ of the worse-anchored of the two arms (see the emit site below).

For an UNSTRANDED library the aligned strand carries no information, so `s_strand` is
taken from a gene annotation instead when one is configured (`--annotation`); see the
"Gene-annotation strand inference" block below.

Coordinates are 0-based, s_start<s_end (circle_start / circle_end convention of truth_bsj.tsv).
Dependency-free (stdlib only). Reads SAM from a file arg or stdin.
"""
import sys, re, os, gzip, argparse

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

# ── Splice-motif strand inference (unstranded libraries only) ────────────────
# The spliceosome leaves its signature in the genome at the two breakpoints, and
# that signature is strand-specific, so it names the junction's strand outright --
# no annotation, and (unlike the host-gene rule below) it is right about ANTISENSE
# circles, which are transcribed opposite to the gene they sit inside.
#
# With `left`/`right` the circle's 0-based start and end, the two dinucleotides are
#     donor    = genome[right : right+2]      intron base pair just past the circle
#     acceptor = genome[left-2 : left]        intron base pair just before it
# On the transcript's own strand a major-spliceosome intron reads GT..AG and a minor
# (U12) intron AT..AC. Read off the PLUS strand of the genome, a minus-strand circle
# therefore shows the reverse complement, with donor and acceptor swapping ends:
#
#     (donor, acceptor)   strand   intron
#       GT   AG              +     U2  (major)
#       AT   AC              +     U12 (minor)
#       CT   AC              -     U2  rc
#       GT   AT              -     U12 rc
#
# All four pairs are distinct, so the PAIR is unambiguous even though GT and AC each
# appear on both strands. Anything else is a non-canonical junction and yields '.'.
MOTIF_STRAND = {("GT", "AG"): "+", ("AT", "AC"): "+",
                ("CT", "AC"): "-", ("GT", "AT"): "-"}


class GenomeSeq:
    """Dinucleotide lookups on a reference FASTA.

    Prefers faidx random access -- a junction needs 4 bases, so reading 3 GB into
    every concurrent detect job to get them would be absurd. Falls back to an
    in-memory parse when there is no .fai or the FASTA is gzipped (a plain gzip is
    not seekable; only bgzf is), which is the small-genome / benchmark case.
    """

    def __init__(self, path):
        self.path = str(path)
        self.idx, self.mem, self.fh = None, None, None
        fai = self.path + ".fai"
        if not self.path.endswith(".gz") and os.path.exists(fai):
            self.idx = {}
            with open(fai) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 5:
                        # name -> (length, offset, bases-per-line, bytes-per-line)
                        self.idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
            self.fh = open(self.path, "rb")
        else:
            self.mem = {}
            opener = gzip.open if self.path.endswith(".gz") else open
            name, buf = None, []
            with opener(self.path, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if name:
                            self.mem[name] = "".join(buf)
                        name, buf = line[1:].split()[0], []
                    else:
                        buf.append(line.strip())
            if name:
                self.mem[name] = "".join(buf)

    def fetch(self, chrom, start, end):
        """0-based half-open slice, uppercased; '' if out of range or unknown contig."""
        if start < 0:
            return ""
        if self.mem is not None:
            seq = self.mem.get(chrom)
            return seq[start:end].upper() if seq and end <= len(seq) else ""
        rec = self.idx.get(chrom)
        if rec is None:
            return ""
        length, offset, linebases, linewidth = rec
        if end > length:
            return ""
        # skip the newline bytes that fall inside the requested span
        self.fh.seek(offset + start // linebases * linewidth + start % linebases)
        raw = self.fh.read((end - start) + (end - start) // linebases + 2)
        return raw.decode("ascii", "replace").replace("\n", "").replace("\r", "")[:end - start].upper()


def motif_strand(genome, chrom, left, right, slop=2):
    """Strand from the splice motifs at the breakpoints, or '.' if non-canonical.

    Microhomology makes a breakpoint ambiguous over a few bases, and the detector
    left-aligns it, so the spliceosomal pair can sit a shift away from the reported
    coordinate. `slop` walks outward from 0 and takes the first shift that reads
    canonical; 0 checks the reported position only.
    """
    for d in range(0, slop + 1):
        for shift in ((0,) if d == 0 else (d, -d)):
            donor = genome.fetch(chrom, right + shift, right + shift + 2)
            acceptor = genome.fetch(chrom, left + shift - 2, left + shift)
            hit = MOTIF_STRAND.get((donor, acceptor))
            if hit:
                return hit
    return "."


# ── Gene-annotation strand inference (unstranded libraries only) ─────────────
# An unstranded library gives no information about which strand a junction came
# from: R1 is sense half the time and antisense the other half, so the ALIGNMENT
# strand of a supporting read is a coin flip. That is not merely uninformative, it
# is actively destructive -- one real junction gets reported at '+' from some of its
# reads and at '-' from the rest, and because the pool dedup keys on the 4-field
# coordinate, it then survives as two half-supported rows.
#
# The fix is to derive the strand from something that does not depend on which mate
# was sequenced. A circRNA is spliced out of a host gene, so the host gene's strand
# IS the junction's strand; a BSJ contained within exactly one annotated gene
# inherits it, and both read orientations then agree. Unresolvable cases return '.'
# rather than a guess -- '.' is still ORIENTATION-INDEPENDENT, so the two halves
# collapse to one row even when the strand itself is unknown.
#
# Genome-scale lookup: genes are bucketed into fixed bins and a junction only probes
# the bin holding its left breakpoint, since containment requires the gene to cover
# it. hg38 has ~60k genes over ~23k bins, so a probe scans a handful of intervals.
GENE_BIN = 1 << 17   # 131072 bp
GENE_ID_RE = re.compile(r'gene_id[ =]"?([^";]+)"?')


def load_gene_strands(path):
    """Parse a GTF/GFF into {chrom: {bin: [(start0, end, strand), ...]}}.

    Uses `gene` features. A GTF that carries none (some GENCODE-derived subsets and
    most GFF3-converted files) falls back to the span of each gene_id's transcript or
    exon lines, so the annotation does not have to be in one particular dialect.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    genes, spans = [], {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[6] not in ("+", "-"):
                continue
            try:
                start, end = int(f[3]) - 1, int(f[4])   # GTF is 1-based inclusive
            except ValueError:
                continue
            if f[2] == "gene":
                genes.append((f[0], start, end, f[6]))
            elif f[2] in ("transcript", "exon"):
                m = GENE_ID_RE.search(f[8])
                if m:
                    key = (f[0], m.group(1), f[6])
                    prev = spans.get(key)
                    spans[key] = (min(prev[0], start), max(prev[1], end)) if prev else (start, end)
    if not genes:
        genes = [(c, lo, hi, st) for (c, _gid, st), (lo, hi) in spans.items()]

    index = {}
    for chrom, start, end, strand in genes:
        per_chrom = index.setdefault(chrom, {})
        for b in range(start // GENE_BIN, (end - 1) // GENE_BIN + 1):
            per_chrom.setdefault(b, []).append((start, end, strand))
    return index, len(genes)


def gene_strand(index, chrom, left, right):
    """Strand of the single annotated gene containing BOTH breakpoints, else '.'.

    Containment, not overlap: a circle is spliced out of its host, so a gene that
    covers only one breakpoint is a neighbour, not the host. Two containing genes on
    opposite strands (a nested or antisense pair) are ambiguous and also give '.'.
    """
    hits = index.get(chrom, {}).get(left // GENE_BIN)
    if not hits:
        return "."
    strands = {st for (lo, hi, st) in hits if lo <= left and right <= hi}
    return strands.pop() if len(strands) == 1 else "."


class StrandInferrer:
    """Junction strand for an UNSTRANDED library: splice motif first, host gene second.

    Motif leads because it is direct physical evidence of the splice reaction at these
    exact breakpoints, it needs no annotation at all, and it is the only one of the two
    that can be right about an ANTISENSE circle -- one transcribed opposite to the gene
    it sits inside reads CT/AC where the host-gene rule can only report the host's '+'.
    Measured on the benchmark's unstranded library, that is worth 6/6 on minus-strand
    truth junctions against 4/6 for the annotation alone. The annotation then covers what
    the motif cannot: non-canonical junctions, which are real but carry no spliceosomal
    signature. Together 0.957 strand accuracy, against 0.872 and 0.851 alone.

    Results are memoised. The strand is a property of the junction, but inference runs
    once per SUPPORTING READ, and a well-covered junction has hundreds -- without the
    cache the same four bases get seeked off disk hundreds of times.
    """

    def __init__(self, genome=None, gene_index=None, motif_slop=2):
        self.genome = genome
        self.gene_index = gene_index
        self.motif_slop = motif_slop
        self._cache = {}
        self.n_motif = self.n_gene = self.n_none = 0

    def strand(self, chrom, left, right):
        key = (chrom, left, right)
        hit = self._cache.get(key)
        if hit is None:
            hit = self._infer(chrom, left, right)
            self._cache[key] = hit
        return hit

    def _infer(self, chrom, left, right):
        if self.genome is not None:
            st = motif_strand(self.genome, chrom, left, right, self.motif_slop)
            if st != ".":
                self.n_motif += 1
                return st
        if self.gene_index is not None:
            st = gene_strand(self.gene_index, chrom, left, right)
            if st != ".":
                self.n_gene += 1
                return st
        self.n_none += 1
        return "."

    def summary(self):
        total = self.n_motif + self.n_gene + self.n_none
        return (f"{total:,} distinct junctions: {self.n_motif:,} by splice motif, "
                f"{self.n_gene:,} by host gene, {self.n_none:,} unresolved ('.')")


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
           inferrer=None,
           strandedness="unstranded",
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
            # `strand` so far is the ALIGNMENT strand of the read. For a reverse-stranded
            # library (dUTP / TruSeq Stranded) R1 is antisense to the transcript, so the
            # junction's real strand is the opposite. R2 is reverse-complemented upstream
            # in rename_read_ids_pe, so both mates arrive in R1's orientation and one flip
            # is correct for all of them. forward: already the transcript strand.
            # unstranded: the aligned strand is a coin flip, so infer it from the
            # genome instead -- splice motif first, host gene second (StrandInferrer).
            # With no inferrer there is nothing better to report than the aligned
            # strand, which is why an unstranded library warns about that.
            if strandedness == "reverse":
                out_strand = "-" if strand == "+" else "+"
            elif strandedness == "unstranded" and inferrer is not None:
                out_strand = inferrer.strand(prim["chrom"], left, right)
            else:
                out_strand = strand
            # MAPQ of the LIMITING arm. A junction is only as well placed as its worse-
            # anchored side: a read whose donor arm is unique (60) but whose acceptor arm
            # is multi-mapping (0) does not localise the breakpoint, so min() is the
            # honest summary. Carried per-read; bsj_summary.py takes the max across the
            # reads supporting a junction (its best evidence) for the confidence tier.
            arm_mapq = min(a["mapq"], b["mapq"])
            yield (f"{qname}\t{bsj_read_pos}\t{prim['chrom']}\t{left}\t{right}\t"
                   f"{out_strand}\t{ovl}\t{arm_mapq}\n")

HEADER = ("q_id\tbsj_coord\tchromosome\ts_start\ts_end\ts_strand\t"
          "overlap_trimmed_nt\tmapq\n")


def write_table(sam_fh, out_fh, min_anchor, min_mapq, max_gap_size,
                max_read_gap, max_read_ovlp, strandedness="unstranded",
                inferrer=None):
    """Stream SAM records in, bsj_cor rows out. Returns the row count."""
    out_fh.write(HEADER)
    n = 0
    for line in detect(sam_fh, min_anchor, min_mapq, max_gap_size, inferrer,
                       strandedness, max_read_gap, max_read_ovlp):
        out_fh.write(line)
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="SA-tag BSJ detector for minimap2 SAM")
    ap.add_argument("--sam", default="-", help="SAM file, or '-' for stdin")
    ap.add_argument("--out", default="-", help="output bsj_cor.tsv, or '-' for stdout")
    ap.add_argument("--min-anchor", type=int, default=30)
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--genome", default=None,
                    help="reference FASTA; for an UNSTRANDED library its splice motifs "
                         "give each junction's strand (faidx random access when a .fai "
                         "sits beside it)")
    ap.add_argument("--annotation", default=None,
                    help="GTF/GFF (.gz ok), the FALLBACK for an unstranded junction whose "
                         "motif is non-canonical: the strand of the single annotated gene "
                         "containing both breakpoints. Unresolved junctions get '.'")
    ap.add_argument("--motif-slop", type=int, default=2,
                    help="bp to shift both breakpoints looking for a canonical motif. "
                         "Microhomology makes a breakpoint ambiguous over a few bases "
                         "and the detector left-aligns it, so the spliceosomal pair can "
                         "sit a shift away (default 2; 0 = reported position only)")
    ap.add_argument("--max-gap-size", type=int, default=100000)
    ap.add_argument("--max-read-gap", type=int, default=15)
    ap.add_argument("--max-read-ovlp", type=int, default=15)
    ap.add_argument("--strandedness", default="unstranded",
                    help="forward|reverse|unstranded, or a detect_strandedness TSV path. "
                         "reverse flips the junction strand (R1 antisense; dUTP/TruSeq)")
    a = ap.parse_args()
    fh = sys.stdin if a.sam == "-" else open(a.sam)
    out = sys.stdout if a.out == "-" else open(a.out, "w")
    n = 0
    try:
        strandedness = resolve_strandedness(a.strandedness)
        inferrer = build_inferrer(a.genome, a.annotation, strandedness,
                                  a.motif_slop, sys.stderr.write)
        n = write_table(fh, out, a.min_anchor, a.min_mapq, a.max_gap_size,
                        a.max_read_gap, a.max_read_ovlp,
                        strandedness, inferrer)
        if inferrer is not None:
            sys.stderr.write(f"[sa_detect] strand inference -- {inferrer.summary()}\n")
    finally:
        if fh is not sys.stdin: fh.close()
        if out is not sys.stdout: out.close()
    sys.stderr.write(f"[sa_detect] {n} BSJ read-support rows\n")


def resolve_strandedness(value):
    """Accept a literal (forward/reverse/unstranded) or a detect_strandedness TSV.

    The verdict has to arrive as an INPUT FILE rather than a param: Snakemake resolves
    params during DAG construction, before detect_strandedness has run, so a param
    could never carry a value another rule computes.
    """
    if value in ("forward", "reverse", "unstranded"):
        return value
    with open(value) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
    verdict = dict(zip(header, row))["verdict"]
    if verdict not in ("forward", "reverse", "unstranded"):
        raise SystemExit(f"unusable strandedness verdict {verdict!r} in {value}")
    return verdict


def build_inferrer(genome_path, annotation, strandedness, motif_slop, log):
    """Assemble the unstranded-library strand inferrer, or None where unused.

    Nothing is loaded for forward/reverse libraries: their junction strand is already
    determined by the read, so the genome and annotation would cost load time and
    change no output. An unstranded library with NEITHER keeps the old behaviour (the
    aligned strand) and warns, because quiet degradation there looks like data rather
    than a missing setting.
    """
    if strandedness != "unstranded":
        if annotation:
            log("[sa_detect] annotation ignored: only used for unstranded libraries\n")
        return None

    genome = gene_index = None
    if genome_path and os.path.exists(str(genome_path)):
        genome = GenomeSeq(genome_path)
        log(f"[sa_detect] splice-motif strand inference from {genome_path} "
            f"({'faidx' if genome.idx is not None else 'in-memory'}"
            f"{f', slop {motif_slop}bp' if motif_slop else ''})\n")
    else:
        log(f"[sa_detect] WARNING: genome {genome_path!r} not readable -- no splice-motif "
            f"strand inference\n")
    if annotation:
        gene_index, n_genes = load_gene_strands(annotation)
        log(f"[sa_detect] host-gene fallback: {n_genes:,} genes from {annotation} "
            f"({len(gene_index)} contigs)\n")

    if genome is None and gene_index is None:
        log("[sa_detect] WARNING: unstranded library with neither a readable genome nor "
            "an `annotation` -- junction strand falls back to the ALIGNED strand, which "
            "is arbitrary for such a library and splits junctions across both strands\n")
        return None
    return StrandInferrer(genome, gene_index, motif_slop)


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
    strandedness = resolve_strandedness(str(smk.input.strandedness))
    annotation = getattr(smk.params, "annotation", None) or None
    genome = getattr(smk.params, "genome", None) or None
    motif_slop = int(getattr(smk.params, "motif_slop", 2) or 0)
    out_path = str(smk.output.tsv)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    with open(str(smk.log[0]), "w") as log_fh:
        cmd = ["samtools", "view", "-@", str(smk.threads), "-F", "0x900",
               str(smk.input.bam)]
        log_fh.write(f"[sa_detect] {' '.join(cmd)}\n")
        log_fh.write(f"[sa_detect] strandedness: {strandedness}"
                     f"{' (junction strand flipped)' if strandedness == 'reverse' else ''}\n")
        inferrer = build_inferrer(genome, annotation, strandedness,
                                  motif_slop, log_fh.write)
        log_fh.flush()
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log_fh, text=True)
        try:
            with gzip.open(out_path, "wt") as out_fh:
                n = write_table(proc.stdout, out_fh,
                                p.min_anchor, p.min_mapq, p.max_gap_size,
                                p.max_read_gap, p.max_read_ovlp, strandedness,
                                inferrer)
        finally:
            proc.stdout.close()
            rc = proc.wait()
        if rc != 0:
            # Without this the rule would "succeed" on a truncated BAM: the pipeline's
            # old `samtools | python` form hid samtools' exit status the same way.
            raise SystemExit(f"samtools view failed (exit {rc}) -- see {smk.log[0]}")
        if inferrer is not None:
            log_fh.write(f"[sa_detect] strand inference -- {inferrer.summary()}\n")
        log_fh.write(f"[sa_detect] {n} BSJ read-support rows\n")


if "snakemake" in globals():
    _from_snakemake(globals()["snakemake"])   # noqa: F821  (injected by script:)
elif __name__ == "__main__":
    main()
