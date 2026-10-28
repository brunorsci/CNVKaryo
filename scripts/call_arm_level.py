#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Whole-chromosome and arm-level aneuploidy calling from per-bin log2 ratios.

This reproduces what a conventional karyogram reports -- whole-chromosome
gains and losses, and gains or losses of a whole arm -- from read depth
alone. It is deliberately NOT segment-based. CBS segmentation optimises
breakpoint placement, so on noisy data it emits many short, high-variance
segments, and calling copy number per segment reports an aneuploidy
wherever the noise happens to cross a threshold. Measured on this
pipeline's own validation panel, the segment path produced 156
non-neutral bands across 17 chromosomes on a constitutional trisomy 21
sample, which has exactly one true abnormality.

Design follows KaryoScan (Maxwell et al. 2017, bioRxiv 204719), the
published method for karyotype detection from exome data, on three
points:

  1. The chromosome is tested as a whole first, and only then its arms.
     A conventional karyotype reports +21, not dup(21q).

  2. The baseline for chromosome i is the median over every OTHER
     autosome. KaryoScan excludes the tested chromosome explicitly "to
     avoid dampening the signal of a trisomy"; including it lets a large
     event pull the very baseline it is measured against.

  3. A call must be significant against the sample's own noise, not
     against a fixed log2 cutoff. Noise is estimated block-to-block
     WITHIN each region and pooled across the genome, so a true
     region-wide shift moves every block together and does not inflate
     the estimate. Estimating it between regions instead collapses on a
     heavily aneuploid genome, where most regions are altered and the
     spread is mostly real signal (measured: on HCC1395, 28 of 39 arms
     altered, a between-arm estimate gave sigma=0.375 and recovered 6 of
     41 documented events; the within-region estimate gave sigma=0.188
     and recovered 22).

KaryoScan additionally restricts bins to 45-55% GC. That was tested here
and rejected: on this pipeline's bins it amplified the chromosome-scale
artefact it was meant to suppress (HCC1395BL chr6 went from -0.142 to
-0.219, chr16 from -0.139 to -0.201, and the genome baseline shifted
from -0.011 to -0.047). CNVkit already applies a GC correction in `fix`,
so a second GC restriction on top of it is not the same operation it is
in the CLAMMS/KaryoScan pipeline.

KaryoScan's other signal -- B-allele balance, where a trisomy drives
heterozygous sites towards 1/3 and 2/3 and a monosomy removes
heterozygosity altogether -- is orthogonal to read depth and is what
lets it separate copy-neutral LOH from true dosage change. It is not
implemented here; this tool is depth-only, and that limit belongs in the
methods.

Usage:
    python3 call_arm_level.py --cnr sample.cnr --cytoband cytoBand.txt \
        --sex-file sample_sex.txt [--declared-sex Female] --sample ID \
        --out-arms arms.tsv --out-karyotype karyotype.txt
"""
import argparse
import csv
import itertools
import math
import os
import statistics
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnvkaryo_plot_common import CHR_ORDER, CHR_LENGTHS, load_cytoband, ploidy_expected

AUTOSOMES = [c for c in CHR_ORDER if c not in ("chrX", "chrY")]

# Short arms of the acrocentric chromosomes: satellite and rDNA repeat,
# no unique exonic sequence, and not scored in conventional cytogenetics
# either. The few bins that map there are repeat artefacts -- on the
# complex-karyotype AML, 21p held 93 bins against 3006 on 21q and called
# a spurious dup(21p), and 13p/14p/15p held 1-4 bins each at log2 -1.5
# to -5.6. Reported for audit, never called.
ACROCENTRIC = {"chr13", "chr14", "chr15", "chr21", "chr22"}

ARMS_HEADER = [
    "sample", "chrom", "region", "start", "end", "n_bins", "n_blocks",
    "median_log2", "expected_level_log2", "delta_log2", "z_score",
    "frac_bins_supporting", "expected_ploidy", "call", "iscn",
]


def parse_sex(sex_file):
    with open(sex_file) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    rec = dict(zip(lines[0].split("\t"), lines[-1].split("\t")))
    return rec.get("sex")


def centromere(chrom, cytobands):
    """(acen_start, acen_end), or None when the cytoband file has no acen
    entry for this contig."""
    acen = [b for b in cytobands.get(chrom, []) if b[3] == "acen"]
    return (min(b[0] for b in acen), max(b[1] for b in acen)) if acen else None


def load_bins(cnr_path, min_weight):
    """chrom -> list of (midpoint, log2).

    Zero-depth bins are dropped: they carry no evidence and their log2 is
    an artefact of dividing by the reference, which is what poisons a
    genome-wide median.

    Weight filtering is off by default. CNVkit's per-bin weight tracks
    depth and GC, so any non-zero cut removes bins non-randomly and
    shifts the baseline rather than just discarding noise. Measured on
    Trisomy21_G3143 (19.6x mean target depth, every bin under weight
    0.06): at 0.1 no bin survives and the sample is uncallable; at 0.01
    the baseline falls from -0.013 to -0.295 and four spurious calls
    appear; at 0 the only call is the documented one. This statistic
    takes a median over thousands of bins, so it does not need per-bin
    reliability -- it needs an unbiased baseline.
    """
    bins = {c: [] for c in CHR_ORDER}
    with open(cnr_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            chrom = r["chromosome"]
            if chrom not in bins or float(r["depth"]) <= 0:
                continue
            if min_weight and float(r["weight"]) < min_weight:
                continue
            bins[chrom].append(((int(r["start"]) + int(r["end"])) / 2, float(r["log2"])))
    return bins


def block_medians(region_bins, block_bp, min_bins_block):
    """Median log2 per fixed-width block, for the within-region noise
    estimate. Blocks are genomic so each carries a comparable amount of
    sequence regardless of exon density."""
    blocks = {}
    for mid, l2 in region_bins:
        blocks.setdefault(int(mid // block_bp), []).append(l2)
    return [statistics.median(v) for v in blocks.values() if len(v) >= min_bins_block]


def iscn_label(chrom, region, call):
    """Conventional shorthand: +8 / -7 for a whole chromosome, del(5q) /
    dup(1q) for a single arm."""
    if call not in ("gain", "loss"):
        return ""
    n = chrom.replace("chr", "")
    if region == "whole":
        return f"{'+' if call == 'gain' else '-'}{n}"
    return f"{'dup' if call == 'gain' else 'del'}({n}{region})"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cnr", required=True)
    ap.add_argument("--cytoband", required=True)
    ap.add_argument("--sex-file", required=True)
    ap.add_argument("--declared-sex", default=None)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-arms", required=True)
    ap.add_argument("--out-karyotype", required=True)
    ap.add_argument("--min-shift", type=float, default=0.25,
                    help="Minimum |median - expected level| in log2. Default 0.25 is set "
                         "from this pipeline's measured negative-control floor, not from "
                         "theory: across the karyotypically normal samples the largest "
                         "whole-chromosome deviation on any autosome was 0.131, while a "
                         "constitutional trisomy 21 measured +0.609 and +0.555 against a "
                         "theoretical +0.585. A single-copy change in a pure diploid sits "
                         "at +0.585 or -1.0, so 0.25 keeps full sensitivity to non-mosaic "
                         "events while staying clear of the floor. Events below it -- "
                         "low-level mosaics, small subclones -- are outside what read "
                         "depth alone resolves here.")
    ap.add_argument("--min-z", type=float, default=4.0,
                    help="Minimum |z| of the shift against the sample's own "
                         "block-to-block noise.")
    ap.add_argument("--min-fraction", type=float, default=0.55,
                    help="Minimum fraction of the region's bins on the same side of the "
                         "expected level -- a handful of extreme bins cannot carry a "
                         "chromosome on their own.")
    ap.add_argument("--min-bins", type=int, default=50,
                    help="Regions with fewer informative bins are reported, never called.")
    ap.add_argument("--min-blocks", type=int, default=2,
                    help="Regions spanning fewer blocks have no usable noise estimate.")
    ap.add_argument("--block-mb", type=float, default=10.0,
                    help="Block width in Mb for the within-region noise estimate.")
    ap.add_argument("--min-bins-block", type=int, default=20,
                    help="Blocks with fewer bins are ignored in the noise estimate.")
    ap.add_argument("--min-weight", type=float, default=0.0,
                    help="Drop bins below this CNVkit weight. Default 0 -- see load_bins().")
    ap.add_argument("--call-sex-chroms", action="store_true",
                    help="Also call chrX/chrY. Off by default. Exome capture on chrY is too "
                         "sparse and repetitive to read this way (median log2 ranged -5.36 "
                         "to +1.19 across samples that should be equivalent), and chrX "
                         "depends on a sex assignment that the coverage itself can "
                         "contradict (HCC1395BL: donor recorded female, chrX observed "
                         "haploid at -0.95). Both need a sex-matched reference and a "
                         "clinically supplied sex before they can be trusted.")
    args = ap.parse_args()

    sex = args.declared_sex or parse_sex(args.sex_file)
    if sex not in ("Male", "Female"):
        sys.exit(f"Invalid sex {sex!r}: must be 'Male' or 'Female'")

    cytobands = load_cytoband(args.cytoband)
    bins = load_bins(args.cnr, args.min_weight)
    block_bp = args.block_mb * 1e6

    per_chrom = {c: [l2 for _, l2 in bins[c]] for c in AUTOSOMES}
    if not any(per_chrom.values()):
        sys.exit(f"[{args.sample}] No informative autosomal bins (depth>0, "
                 f"weight>={args.min_weight}) -- nothing to call against.")

    targets = [c for c in CHR_ORDER
               if c in AUTOSOMES or (args.call_sex_chroms and c in ("chrX", "chrY"))]

    candidates = []
    for chrom in targets:
        end = CHR_LENGTHS[chrom]
        cen = centromere(chrom, cytobands)
        candidates.append((chrom, "whole", 0, end))
        if cen:
            if chrom not in ACROCENTRIC:
                candidates.append((chrom, "p", 0, cen[0]))
            else:
                candidates.append((chrom, "p", 0, cen[0]))
            candidates.append((chrom, "q", cen[1], end))
        else:
            candidates.append((chrom, "q", 0, end))

    rows = []
    for chrom, region, start, end in candidates:
        expected = ploidy_expected(chrom, sex)
        region_bins = [(mid, l2) for mid, l2 in bins[chrom] if start <= mid < end]
        n = len(region_bins)
        # Leave-one-chromosome-out: the tested chromosome is excluded so a
        # whole-chromosome event cannot dampen its own signal (KaryoScan).
        others = list(itertools.chain.from_iterable(
            v for c, v in per_chrom.items() if c != chrom))
        baseline = statistics.median(others)
        # A sex chromosome at its own expected ploidy sits a fixed offset
        # from the autosomal baseline; comparing it directly to that
        # baseline calls -X on every male sample.
        level = baseline + math.log2(expected / 2) if expected else baseline
        row = dict(chrom=chrom, region=region, start=start, end=end, n_bins=n,
                   n_blocks=0, median_log2="", expected_level_log2=round(level, 4),
                   delta_log2="", z_score="", frac_bins_supporting="",
                   expected_ploidy=expected, call="no_data", iscn="")
        if chrom in ACROCENTRIC and region == "p":
            row["call"] = "uninformative"
            if n:
                row["median_log2"] = round(statistics.median([l2 for _, l2 in region_bins]), 4)
        elif n:
            med = statistics.median([l2 for _, l2 in region_bins])
            delta = med - level
            d = 1 if delta > 0 else -1
            bmed = block_medians(region_bins, block_bp, args.min_bins_block)
            row.update(median_log2=round(med, 4), delta_log2=round(delta, 4),
                       n_blocks=len(bmed), call="neutral",
                       frac_bins_supporting=round(
                           sum(1 for _, l2 in region_bins if (l2 - level) * d > 0) / n, 3))
            row["_bmed"] = bmed
            row["_delta"] = delta
        rows.append(row)

    # One per-sample technical sigma, pooled from every arm's own block
    # deviations. Pooling keeps it stable on arms with few blocks;
    # measuring within arms keeps real signal out of it. Whole-chromosome
    # rows are excluded from the pool because a chromosome with one
    # altered arm has genuine block-to-block structure.
    devs = []
    for row in rows:
        bmed = row.get("_bmed") or []
        if len(bmed) >= 2 and row["region"] != "whole":
            c = statistics.median(bmed)
            devs.extend(abs(x - c) for x in bmed)
    sigma = (1.4826 * statistics.median(devs)) if devs else 0.0
    sigma = sigma or 1e-9

    for row in rows:
        if row["call"] != "neutral" or row["n_bins"] < args.min_bins:
            continue
        if row["n_blocks"] < args.min_blocks:
            continue
        z = row["_delta"] / (sigma / math.sqrt(row["n_blocks"]))
        row["z_score"] = round(z, 2)
        if (abs(z) >= args.min_z and abs(row["_delta"]) >= args.min_shift
                and row["frac_bins_supporting"] >= args.min_fraction):
            row["call"] = "gain" if row["_delta"] > 0 else "loss"

    # A whole-chromosome call is only valid when the chromosome behaves as
    # a whole. Testing the whole chromosome as one region lets a mixture
    # of one neutral arm and one altered arm cross the threshold: on
    # BV-173, chr9p sits at -0.06 and chr9q at -0.85, and their combined
    # median of -0.61 was called -9, which then subsumed and erased the
    # correct del(9q). A whole-chromosome call is therefore kept only
    # when every callable arm agrees with it, and demoted otherwise so
    # the arm calls stand. This is what a conventional karyotype does:
    # -9 means the chromosome is gone, del(9q) means an arm is.
    by_chrom = {}
    for row in rows:
        by_chrom.setdefault(row["chrom"], []).append(row)
    for chrom, group in by_chrom.items():
        whole = next((r for r in group if r["region"] == "whole"), None)
        if not whole or whole["call"] not in ("gain", "loss"):
            continue
        arms = [r for r in group
                if r["region"] in ("p", "q") and r["call"] != "uninformative"
                and r["n_bins"] >= args.min_bins]
        agree = arms and all(r["call"] == whole["call"] for r in arms)
        if agree:
            for r in arms:
                r["call"] = "subsumed"
        else:
            whole["call"] = "split_by_arm"

    for row in rows:
        row["iscn"] = iscn_label(row["chrom"], row["region"], row["call"])
        row.pop("_bmed", None)
        row.pop("_delta", None)
        row["iscn"] = iscn_label(row["chrom"], row["region"], row["call"])
        row.pop("_bmed", None)
        row.pop("_delta", None)

    with open(args.out_arms, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ARMS_HEADER, delimiter="\t")
        w.writeheader()
        for row in rows:
            row["sample"] = args.sample
            w.writerow(row)

    findings = []
    for row in rows:
        if row["iscn"] and row["iscn"] not in findings:
            findings.append(row["iscn"])
    called_bp = sum(row["end"] - row["start"] for row in rows
                    if row["call"] in ("gain", "loss"))
    genome_bp = sum(CHR_LENGTHS[c] for c in AUTOSOMES)
    with open(args.out_karyotype, "w") as f:
        f.write(f"Sample: {args.sample}\n")
        f.write(f"Sex used: {sex}\n")
        f.write(f"Technical noise sigma (block-to-block log2): {sigma:.4f}\n")
        f.write(f"Thresholds: |delta|>={args.min_shift} |z|>={args.min_z} "
                f"frac>={args.min_fraction} bins>={args.min_bins}\n")
        f.write(f"Fraction of autosomal sequence in a called region: "
                f"{called_bp / genome_bp:.4f}\n")
        f.write(f"Findings: {', '.join(findings) if findings else 'None'}\n")

    print(f"[{args.sample}] sex={sex} sigma={sigma:.4f} findings={findings or 'None'}")


if __name__ == "__main__":
    main()
