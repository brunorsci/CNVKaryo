#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-sample noise, measured the way a reader judges a copy-number figure.

Published exome and array CNV figures are legible because the trace sits
in a tight band around zero and an event displaces it visibly. The
quantity that decides legibility is therefore not per-bin scatter -- exon
capture makes that large and it cannot be fixed -- but the spread of
medians over fixed genomic windows, which is what the eye integrates when
it follows a trace.

Reported as the MAD of window medians about the autosomal baseline. A
single-copy gain in a pure diploid displaces the trace by 0.585 log2, so
this number states directly how visible such an event will be: at MAD
0.05 it stands 12 window-noise units clear, at MAD 0.28 only 2.

Windows with too few bins are dropped: an exon-sparse window's median is
not comparable to a dense one's and would inflate the estimate.

Usage:
    python3 window_noise_qc.py --cnr a.cnr [b.cnr ...] --out qc.tsv
"""
import argparse
import csv
import os
import statistics
import sys

HEADER = ["sample", "window_mb", "n_windows", "baseline_log2", "noise_mad",
          "visibility_single_copy_gain", "n_bins_used", "frac_bins_zero_depth"]


def window_medians(cnr_path, window_bp, min_bins, exclude):
    rows_total = used = zero = 0
    windows = {}
    with open(cnr_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows_total += 1
            if float(r["depth"]) <= 0:
                zero += 1
                continue
            chrom = r["chromosome"]
            if chrom in exclude:
                continue
            mid = (int(r["start"]) + int(r["end"])) / 2
            windows.setdefault((chrom, int(mid // window_bp)), []).append(float(r["log2"]))
            used += 1
    medians = {k: statistics.median(v) for k, v in windows.items() if len(v) >= min_bins}
    return medians, used, (zero / rows_total if rows_total else 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cnr", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window-mb", type=float, default=1.0)
    ap.add_argument("--min-bins-window", type=int, default=10)
    ap.add_argument("--max-noise", type=float, default=0.05,
                    help="Gate. Above this a sample's trace is too noisy for a "
                         "single-copy event to read off the figure, and the sample "
                         "should be reported as failing QC rather than called.")
    args = ap.parse_args()

    exclude = {"chrX", "chrY"}
    rows = []
    for path in args.cnr:
        sample = os.path.basename(path).replace("_markduplicates.cnr", "").replace(".cnr", "")
        med, used, fz = window_medians(path, args.window_mb * 1e6,
                                       args.min_bins_window, exclude)
        if len(med) < 20:
            print(f"[{sample}] too few usable windows ({len(med)})", file=sys.stderr)
            continue
        baseline = statistics.median(med.values())
        dev = [abs(v - baseline) for v in med.values()]
        mad = 1.4826 * statistics.median(dev)
        rows.append(dict(sample=sample, window_mb=args.window_mb, n_windows=len(med),
                         baseline_log2=round(baseline, 4), noise_mad=round(mad, 4),
                         visibility_single_copy_gain=round(0.585 / mad, 1) if mad else "",
                         n_bins_used=used, frac_bins_zero_depth=round(fz, 4)))

    rows.sort(key=lambda r: r["noise_mad"])
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"{'sample':<34}{'noise(MAD)':>11}{'visibility':>12}{'zero-depth':>12}  gate")
    for r in rows:
        ok = "PASS" if r["noise_mad"] <= args.max_noise else "FAIL"
        print(f"{r['sample'][:33]:<34}{r['noise_mad']:>11.4f}"
              f"{r['visibility_single_copy_gain']:>12}"
              f"{r['frac_bins_zero_depth']:>11.1%}  {ok}")
    n_pass = sum(1 for r in rows if r["noise_mad"] <= args.max_noise)
    print(f"\n{n_pass}/{len(rows)} at or below MAD {args.max_noise} -> {args.out}")


if __name__ == "__main__":
    main()
