#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Detect homozygous deletions as runs of intervals with no reads, where the
reference panel has reads.

A biallelic loss leaves no template, so coverage goes to zero rather than
to half. That is a qualitatively different observation from the log2
shift every other caller in this pipeline works from, and it is thrown
away by them: the log2 of zero coverage is meaningless (CNVkit emits
about -24), and the arm-level caller discards zero-depth intervals by
design because elsewhere they are artefacts of dividing by the
reference.

Measured on the LL-100 panel: GRANTA-519, whose DSMZ entry records
"deletion of p15/p16", has 34 consecutive intervals at exactly 0.0x
across CDKN2A/CDKN2B at 9p21, while sixteen of twenty panel members
carry 27-53x over the same intervals. U-2940, recorded as carrying a
"homozygous microdeletion of SOCS1", shows the same at 16p13. The signal
is unambiguous; nothing was detecting it.

The discrimination that makes this specific is the panel. Zero coverage
alone means either a biallelic deletion or a region the capture never
reaches. Requiring the panel to have solid coverage over the same
intervals separates the two, using the same data the reference is built
from.

Events smaller than --min-intervals are not reported: a single empty
interval is more often a mapping artefact than a deletion.

Usage:
    python3 call_homozygous_deletions.py --coverage sample.targetcoverage.cnn \
        --panel-coverage p1.cnn p2.cnn ... --sample ID --out calls.tsv
"""
import argparse
import csv
import os
import statistics
import sys

HEADER = ["sample", "event_type", "chrom", "start", "end", "n_intervals", "span_bp",
          "sample_depth_median", "panel_depth_median",
          "panel_fraction_covered", "log2_median", "genes", "flag"]

# Regions where absent coverage is expected and is not a copy-number
# finding. They are flagged and kept in the output rather than dropped,
# so the call can be audited, but they are not reported as deletions.
#
# Immunoglobulin and T-cell receptor loci undergo V(D)J recombination,
# which physically excises segments in any lymphoid sample. Every B-cell
# line in this cohort shows IGH at 14q32 and IGL at 22q11 as empty; that
# is normal lymphocyte biology, not a deletion. GRCh38 coordinates.
IG_TCR_LOCI = [
    ("chr14", 105550000, 106900000, "IGH"),
    ("chr22", 22020000, 22930000, "IGL"),
    ("chr2", 88800000, 90300000, "IGK"),
    ("chr14", 21620000, 22560000, "TRA/TRD"),
    ("chr7", 142290000, 142820000, "TRB"),
    ("chr7", 38240000, 38380000, "TRG"),
]

# Genes whose biallelic deletion is a common germline polymorphism
# rather than a finding. GSTM1 null occurs in roughly half of the
# population and UGT2B17 null in a quarter or more.
COMMON_CNV_GENES = {"GSTM1", "GSTT1", "UGT2B17", "LILRA6", "LILRA3",
                    "RHCE", "RHD", "FAM66A", "AMY1A", "NAIP"}


def read_coverage(path):
    """(chrom, start, end) -> (depth, gene), in file order."""
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            try:
                d = float(r["depth"])
            except (TypeError, ValueError, KeyError):
                continue
            rows.append(((r["chromosome"], int(r["start"]), int(r["end"])),
                         d, r.get("gene", "")))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--coverage", required=True,
                    help="This sample's raw target coverage (.targetcoverage.cnn). Raw "
                         "depth, not log2: the ratio is undefined where coverage is zero.")
    ap.add_argument("--panel-coverage", nargs="+", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-sample-depth", type=float, default=0.5,
                    help="Depth at or below which the sample counts as empty. Not exactly "
                         "zero, because a stray mismapped read should not rescue an "
                         "interval that is otherwise absent.")
    ap.add_argument("--min-panel-depth", type=float, default=10.0,
                    help="Median panel depth required over the interval. Below this the "
                         "panel cannot say whether the capture reaches the region, so "
                         "absence in the sample is uninformative.")
    ap.add_argument("--min-panel-fraction", type=float, default=0.5,
                    help="Fraction of panel members that must have reads over the "
                         "interval. A recurrent deletion, such as CDKN2A in lymphoid "
                         "lines, can leave several panel members empty too: four of "
                         "twenty here. The majority must still be covered.")
    ap.add_argument("--max-gap-bp", type=int, default=1000000,
                    help="Intervals further apart than this start a new event. Without "
                         "it, empty intervals at opposite ends of a chromosome merge "
                         "into one span of over a hundred megabases, which is not an "
                         "event but a bookkeeping error.")
    ap.add_argument("--include-sex-chroms", action="store_true",
                    help="Score chrX and chrY. Off by default: a male sample has one X "
                         "against a panel of mixed sex, and chrY is absent from half the "
                         "panel, so coverage differences there are sex, not deletion.")
    ap.add_argument("--ratios", default=None,
                    help="This sample's .cnr. Enables focal amplification calling, which "
                         "needs the normalised ratio rather than raw depth.")
    ap.add_argument("--min-amp-log2", type=float, default=1.5,
                    help="Log2 ratio above which an interval counts as amplified. 1.5 is "
                         "about five copies against two and is far above anything "
                         "whole-chromosome: a trisomy is 0.585. Measured here, the BCL2 "
                         "locus reaches +3.2 in RI-1 and +3.3 in U-2932, both recorded in "
                         "the LL-100 table as 'BCL2 amp', against +0.14 in a line without "
                         "it.")
    ap.add_argument("--min-intervals", type=int, default=3,
                    help="Consecutive empty intervals required. One is more often a "
                         "mapping artefact than a deletion.")
    args = ap.parse_args()

    sample = read_coverage(args.coverage)
    if not sample:
        sys.exit(f"[{args.sample}] no usable intervals in {args.coverage}")

    panel = {}
    for p in args.panel_coverage:
        if os.path.abspath(p) == os.path.abspath(args.coverage):
            continue  # a sample cannot be its own control
        for key, depth, _ in read_coverage(p):
            panel.setdefault(key, []).append(depth)
    if not panel:
        sys.exit(f"[{args.sample}] panel is empty after excluding the sample itself")

    def region_flag(chrom, start, end, genes):
        for c, a, b, name in IG_TCR_LOCI:
            if chrom == c and start < b and end > a:
                return f"{name}_locus_VDJ_recombination"
        if genes and all(g in COMMON_CNV_GENES for g in genes):
            return "common_germline_CNV"
        return ""

    flagged = []
    for key, depth, gene in sample:
        if not args.include_sex_chroms and key[0] in ("chrX", "chrY"):
            continue
        ref = panel.get(key)
        if not ref:
            continue
        covered = sum(1 for d in ref if d > 0) / len(ref)
        if (depth <= args.max_sample_depth
                and statistics.median(ref) >= args.min_panel_depth
                and covered >= args.min_panel_fraction):
            flagged.append((key, depth, statistics.median(ref), covered, gene))

    # Focal amplification, read from the normalised ratio rather than from
    # raw depth: unlike a biallelic deletion, a high-level gain is not a
    # qualitative change in coverage and the reference must be divided
    # out before it can be judged.
    amp_flagged = []
    if args.ratios:
        with open(args.ratios) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                chrom = r["chromosome"]
                if not args.include_sex_chroms and chrom in ("chrX", "chrY"):
                    continue
                try:
                    d, l2 = float(r["depth"]), float(r["log2"])
                except (TypeError, ValueError, KeyError):
                    continue
                if d <= 0 or l2 < args.min_amp_log2:
                    continue
                key = (chrom, int(r["start"]), int(r["end"]))
                amp_flagged.append((key, l2, 0.0, 1.0, r.get("gene", "")))

    # Merge consecutive flagged intervals on the same chromosome.
    def merge(items):
        out, cur = [], []
        for item in items:
            same = cur and item[0][0] == cur[-1][0][0]
            adj = same and (item[0][1] - cur[-1][0][2]) <= args.max_gap_bp
            if adj:
                cur.append(item)
            else:
                if cur:
                    out.append(cur)
                cur = [item]
        if cur:
            out.append(cur)
        return out

    runs, cur = [], []
    for item in flagged:
        same_chrom = cur and item[0][0] == cur[-1][0][0]
        adjacent = same_chrom and (item[0][1] - cur[-1][0][2]) <= args.max_gap_bp
        if adjacent:
            cur.append(item)
        else:
            if cur:
                runs.append(cur)
            cur = [item]
    if cur:
        runs.append(cur)

    rows = []
    for run in runs + [("AMP", r) for r in merge(amp_flagged)]:
        amp = isinstance(run, tuple) and run[0] == "AMP"
        if amp:
            run = run[1]
        if len(run) < args.min_intervals:
            continue
        chrom = run[0][0][0]
        start, end = run[0][0][1], run[-1][0][2]
        genes = sorted({g for *_, g in run if g and g not in ("-", ".")})
        flag = region_flag(chrom, start, end, genes)
        rows.append(dict(
            flag=flag,
            event_type="focal_amplification" if amp else "homozygous_deletion",
            log2_median=round(statistics.median([d for _, d, _, _, _ in run]), 3) if amp else "",
            sample=args.sample, chrom=chrom, start=start, end=end,
            n_intervals=len(run), span_bp=end - start,
            sample_depth_median="" if amp else round(statistics.median([d for _, d, _, _, _ in run]), 3),
            panel_depth_median="" if amp else round(statistics.median([m for _, _, m, _, _ in run]), 1),
            panel_fraction_covered="" if amp else round(statistics.median([c for *_, c, _ in run]), 3),
            genes=",".join(genes[:20])))

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    rep = [r for r in rows if not r["flag"]]
    nd = sum(1 for r in rep if r["event_type"] == "homozygous_deletion")
    na = len(rep) - nd
    print(f"[{args.sample}] {nd} homozygous deletion(s), {na} focal amplification(s), "
          f"{len(rows) - len(rep)} flagged as expected -> {args.out}")
    for r in rep[:14]:
        extra = (f"log2 {r['log2_median']}" if r["event_type"] == "focal_amplification"
                 else f"panel {r['panel_depth_median']}x")
        print(f"    {r['event_type'][:4]} {r['chrom']}:{r['start']}-{r['end']}  "
              f"{r['n_intervals']} int  {extra}  {r['genes'][:52]}")


if __name__ == "__main__":
    main()
