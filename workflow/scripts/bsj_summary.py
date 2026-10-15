#!/usr/bin/env python3
"""
bsj_summary.py — single summary table of all BSJs found across all samples.

Takes the final unique set of polished BSJs (representatives from the cross-sample
pool dedup) and inner-joins each representative back onto the per-sample bsj_cor
coordinate rows, so only junctions that SURVIVED pool deduplication appear.
The sample of origin is dropped (per design): the same junction called in several
samples collapses to one row. `n_samples` = number of distinct samples contributing
a member to that junction's cluster(s), from the pool cluster TSV.

`n_loci` is the coordinate-concordance check. mmseqs clusters the +/-flank junction
SEQUENCES, so paralogs, segmental duplications and shared-arm (alternative) back-
splices can land in one cluster despite sitting at different loci. n_loci counts the
distinct loci among that cluster's members (breakpoints within `locus_tol` bp of each
other, on the same chromosome and strand, count as one). n_loci == 1 is a clean
cluster; n_loci > 1 means the sequence clustering merged junctions that coordinates
say are distinct, so the row's identity is unreliable.

On GBM ~820 of 763k clusters (0.1%) are discordant. Tightening mmseqs_min_seq_id or
mmseqs_min_cov does NOT fix this: measured across a parameter sweep, every setting
destroyed 2-20x more legitimate breakpoint-wobble merges than false merges it removed,
because the two populations are not separable in sequence space. Coordinates are the
discriminator, hence this check.

`n_samples` counts only members CONCORDANT with the row's own coordinate, so a merged
cluster no longer credits one locus with reads from another. For n_loci == 1 rows this
is identical to counting every member.

`n_reads` is how many cluster members (i.e. supporting reads) are concordant with
the row's own coordinate -- the same members n_samples is counted over. `min_reads`
drops rows below a threshold.

Why it exists: a junction seen once is weak evidence, and single-read calls are also
what populate the residual +/- strand duplicates. Measured on the Tasmanian devil
set, 90.9% of coordinates appearing on both strands involve no orphan reads, and
86.8% of THOSE have exactly one read on the minority strand against a median of 13
on the other. Those are not two junctions; they are one shallow observation reported
twice. The effect is worst where detection yield is low: 75 bp devil reads yield
0.27% junctions per SA record against 7.17% for 150 bp GBM reads, because a 75 bp
read must carry the junction between positions 20 and 55 for both arms to clear
min_align_len=20.

Default is 1 (report everything), so enabling the filter is an explicit per-dataset
choice rather than a silent change to the catalogue.

A row is keyed on the 4-field coordinate (chromosome, s_start, s_end, s_strand)
ONLY. `overlap_trimmed_nt` is a property of the supporting READ -- how much
microhomology sa_detect.py consumed to place that read's breakpoint -- not of the
junction, so reads supporting one junction can disagree on it. The column reports
the MAXIMUM value seen, i.e. the widest breakpoint-ambiguity window observed for
that junction.

Representative header format (pool fasta):
    {sample}.{read_id}_BSJ_{bsj_coord}_region_{a}-{b}
so the join key is (q_id = {sample}.{read_id}, bsj_coord).

Stdlib only (matches sa_detect.py — no conda env needed).
"""

import gzip
import re
import argparse
import collections
import logging

# {q_id}_BSJ_{coord}_region_{a}-{b}  (q_id itself contains dots, e.g. ERR10783600.1.18629)
HEADER_RE = re.compile(r"^(?P<qid>.+?)_BSJ_(?P<coord>\d+)_region_\d+-\d+$")


def chrom_key(chrom):
    """Natural sort: chr1 < chr2 < chr10, non-numeric contigs last (alphabetical)."""
    m = re.match(r"^(?:chr)?(\d+)$", chrom)
    if m:
        return (0, int(m.group(1)), "")
    return (1, 0, chrom)


def sample_of(qid):
    return qid.split(".", 1)[0]


def as_int(value):
    """overlap_trimmed_nt as int; a non-numeric field sorts as 0 rather than crashing."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def same_locus(a, b, tol):
    """Two member coordinates describing the same junction, allowing breakpoint slop."""
    return (a[0] == b[0] and a[3] == b[3]
            and abs(a[1] - b[1]) <= tol and abs(a[2] - b[2]) <= tol)


def group_loci(coords, tol):
    """Single-linkage grouping of member coordinates into distinct loci.

    Clusters are tiny (mean ~1.6 members), so the quadratic scan is irrelevant and
    single-linkage chaining cannot run away.
    """
    groups = []
    for c in sorted(coords):
        for g in groups:
            if any(same_locus(c, o, tol) for o in g):
                g.append(c)
                break
        else:
            groups.append([c])
    return groups


def main(cor, pool_fasta, pool_cluster, out, locus_tol=10, min_reads=1):
    """cor: list of per-sample bsj_cor.tsv.gz; the rest are single paths."""

    # ── 1) Index bsj_cor: (q_id, bsj_coord) -> (chrom, s_start, s_end, s_strand, overlap) ──
    idx = {}
    n_cor_rows = 0
    for path in cor:
        with gzip.open(path, "rt") as fh:
            next(fh, None)  # header line
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 7:
                    continue
                idx[(p[0], p[1])] = (p[2], p[3], p[4], p[5], p[6])
                n_cor_rows += 1
    logging.info(f"bsj_cor rows indexed: {n_cor_rows:,} from {len(cor)} sample files")

    # ── 2) Cluster membership: representative -> set of member samples ─────────
    # `mmseqs createtsv` emits REPRESENTATIVE <TAB> MEMBER, so col1 keys the cluster
    # and col2 is the sequence in it. Keying this the other way round (support[t] +=
    # sample_of(q)) looks plausible but is silently useless: a representative is only
    # ever a member of its own cluster, so every lookup returns just its own sample
    # and all cross-sample support is dropped. Sanity check: len(support) must equal
    # the number of records in the pool fasta.
    # rep -> [(member coordinate or None, member sample), ...]. The coordinate is what
    # makes the concordance check possible; None means the member's header did not
    # resolve in bsj_cor, in which case it is counted as concordant rather than dropped.
    support = collections.defaultdict(list)
    n_members = n_unresolved = 0
    with open(pool_cluster) as fh:
        for line in fh:
            q, _, t = line.rstrip("\n").partition("\t")
            if q == "query":  # header
                continue
            if not t:
                continue
            m = HEADER_RE.match(t)
            row = idx.get((m.group("qid"), m.group("coord"))) if m else None
            coord = (row[0], as_int(row[1]), as_int(row[2]), row[3]) if row else None
            if coord is None:
                n_unresolved += 1
            support[q].append((coord, sample_of(t)))
            n_members += 1
    logging.info(f"cluster members read: {n_members:,} (clusters: {len(support):,}"
                 f"{f', unresolved coords: {n_unresolved:,}' if n_unresolved else ''})")

    # ── 3) Inner join representatives -> coordinate rows ───────────────────────
    # rows:     junction coordinate -> set of pool representatives mapping to it
    # overlaps: junction coordinate -> every overlap_trimmed_nt seen for it
    #
    # Keyed on the 4-field coordinate only (see module docstring). Including
    # overlap_trimmed_nt in the key split a junction into one row per distinct
    # overlap its reads happened to report, and divided n_samples between them.
    rows = collections.defaultdict(set)
    overlaps = collections.defaultdict(set)
    n_reps = n_joined = n_missing = 0
    with open(pool_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            h = line[1:].strip()
            n_reps += 1
            m = HEADER_RE.match(h)
            if m is None:
                n_missing += 1
                logging.warning(f"Unparseable representative header: {h}")
                continue
            row = idx.get((m.group("qid"), m.group("coord")))
            if row is None:
                n_missing += 1
                logging.warning(f"Representative not found in bsj_cor index: {h}")
                continue
            n_joined += 1
            coord = row[:4]
            rows[coord].add(h)
            overlaps[coord].add(as_int(row[4]))

    n_disputed = sum(1 for v in overlaps.values() if len(v) > 1)
    logging.info(f"representatives: {n_reps:,} | joined: {n_joined:,} | missing: {n_missing:,}")
    logging.info(f"unique junctions after dropping sample: {len(rows):,}")
    logging.info(f"junctions whose reads disagree on overlap_trimmed_nt: {n_disputed:,} "
                 f"(reporting the max for each)")

    # ── 4) Write summary table ──────────────────────────────────────────────────
    tol = locus_tol
    n_discordant = 0
    n_filtered = 0
    with open(out, "w") as fo:
        fo.write("chromosome\ts_start\ts_end\ts_strand\toverlap_trimmed_nt\t"
                 "n_reads\tn_samples\tn_loci\n")
        for row in sorted(rows, key=lambda r: (chrom_key(r[0]), int(r[1]), int(r[2]))):
            this = (row[0], as_int(row[1]), as_int(row[2]), row[3])
            samples = set()
            member_coords = set()
            n_reads = 0
            for rep in rows[row]:
                for coord, sample in support.get(rep, ()):
                    if coord is None:
                        samples.add(sample)      # unresolvable: keep, do not penalise
                        n_reads += 1
                        continue
                    member_coords.add(coord)
                    if same_locus(coord, this, tol):
                        samples.add(sample)
                        n_reads += 1             # counted over the SAME members as n_samples
            if n_reads < min_reads:
                n_filtered += 1
                continue
            n_loci = len(group_loci(member_coords, tol)) if member_coords else 1
            if n_loci > 1:
                n_discordant += 1
            ovl = max(overlaps[row])
            fo.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{ovl}\t"
                     f"{n_reads}\t{len(samples)}\t{n_loci}\n")

    if min_reads > 1:
        logging.info(f"rows dropped by min_reads={min_reads}: {n_filtered:,} "
                     f"({100*n_filtered/max(len(rows),1):.2f}% of {len(rows):,})")
    logging.info(f"rows whose cluster spans >1 locus (n_loci>1, tol={tol}bp): {n_discordant:,} "
                 f"({100*n_discordant/max(len(rows),1):.2f}%) -- their n_samples counts only "
                 f"concordant members")
    logging.info(f"wrote summary table -> {out}")


def _from_cli():
    """Standalone invocation; kept so the script can be re-run by hand on an existing
    results tree without going through Snakemake."""
    ap = argparse.ArgumentParser(description="Cross-sample BSJ summary table")
    ap.add_argument("--cor", action="append", required=True,
                    help="per-sample bsj_cor.tsv.gz (repeatable)")
    ap.add_argument("--pool-fasta", required=True, help="pool dedup representative fasta")
    ap.add_argument("--pool-cluster", required=True, help="pool dedup cluster tsv")
    ap.add_argument("--out", required=True, help="output summary tsv")
    ap.add_argument("--min-reads", type=int, default=1,
                    help="drop junctions with fewer than this many concordant supporting "
                         "reads (default 1 = no filtering)")
    ap.add_argument("--locus-tol", type=int, default=10,
                    help="breakpoint slop (bp) within which two member coordinates are "
                         "the same locus (default 10; microhomology wobble is typically "
                         "<=10bp, real false merges are far larger)")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main(a.cor, a.pool_fasta, a.pool_cluster, a.out, a.locus_tol, a.min_reads)


# Snakemake's `script:` directive injects a `snakemake` object into the module globals
# and runs the file directly -- there is no argv and __name__ is not reliable across
# versions, so probe for the object rather than testing __name__.
if "snakemake" in globals():
    snakemake = globals()["snakemake"]          # noqa: F821  (injected)
    # `script:` does NOT redirect stdout/stderr to log: the way `shell:` does, so wire
    # logging to the declared log file explicitly or the run is silent.
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s",
                        filename=str(snakemake.log[0]), filemode="w")
    main(
        cor=list(snakemake.input.cor),
        pool_fasta=str(snakemake.input.pool_fasta),
        pool_cluster=str(snakemake.input.pool_cluster),
        out=str(snakemake.output.tsv),
        locus_tol=int(snakemake.params.locus_tol),
        min_reads=int(snakemake.params.min_reads),
    )
elif __name__ == "__main__":
    _from_cli()
