#!/usr/bin/env python3
"""
Builds a kit-agnostic exome target BED from a UCSC refFlat file: the
genome-wide union of all RefSeq exons, merged where they overlap.

Used as the target definition for samples whose actual hybrid-capture kit
is unknown (public data from an external study) -- a real capture kit's
own bait design always differs somewhat from plain exon boundaries, so
this is less precise than a kit-matched BED, but it does not privilege any
one vendor's probe geometry the way reusing another study's kit-specific
BED would (that mismatch is what produced biologically impossible cn
values -- cn=28, cn=67 -- when an unrelated kit's reference was reused
across studies).

Usage:
    python3 refflat_to_exon_bed.py --refflat refFlat.txt --output exons.bed
"""

import argparse
from collections import defaultdict


def parse_exons(refflat_path):
    exons_by_chrom = defaultdict(list)
    with open(refflat_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            chrom = fields[2]
            starts = fields[9].rstrip(",").split(",")
            ends = fields[10].rstrip(",").split(",")
            for s, e in zip(starts, ends):
                if not s or not e:
                    continue
                exons_by_chrom[chrom].append((int(s), int(e)))
    return exons_by_chrom


def merge_intervals(intervals):
    intervals.sort()
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refflat", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    exons_by_chrom = parse_exons(args.refflat)

    n_exons = 0
    with open(args.output, "w") as out:
        for chrom in CHR_ORDER:
            if chrom not in exons_by_chrom:
                continue
            merged = merge_intervals(exons_by_chrom[chrom])
            n_exons += len(merged)
            for start, end in merged:
                out.write(f"{chrom}\t{start}\t{end}\n")

    print(f"Merged exon regions written: {n_exons} (across {len(CHR_ORDER)} standard chromosomes)")
