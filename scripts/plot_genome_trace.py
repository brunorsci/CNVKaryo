#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genome-wide copy-number trace, in the layout used by published array-CGH
and WGS copy-number figures.

One continuous axis across the genome, chromosomes delimited by vertical
rules, log2 on y with a fixed range, and the data drawn as a dense point
cloud rather than as 24 separate panels. A gain or a loss is visible
because the cloud moves off zero, which is the only thing that makes such
a figure evidence rather than an illustration of a call already made.

Two layers, and the distinction is the point of the figure:
  - every informative bin, as faint points, so the reader sees the real
    dispersion of the data and can judge the call;
  - the median over fixed genomic windows, as a solid line, which is what
    the eye integrates when following a published trace. Exon capture
    makes per-bin dispersion large and no normalisation removes it; the
    window median is where the signal lives.

Nothing on the figure is derived from a call. No shading, no arrows, no
ISCN labels. Earlier versions drew the caller's verdict onto the image,
which made the figure circular as evidence: it displayed a conclusion
instead of letting the reader reach one. Calls live in the accompanying
tables; the figure carries the data and the axis.

Reads either CNVkit .cnr or a GATK denoised-copy-ratio table, so the
same figure can be drawn from either caller's output and the two can be
compared on equal terms.

Usage:
    python3 plot_genome_trace.py --ratios sample.cnr --sample ID \
        --out trace.png [--window-mb 1]
"""
import argparse
import csv
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnvkaryo_plot_common import CHR_ORDER, CHR_LENGTHS, load_cytoband
from call_arm_level import ACROCENTRIC, centromere


def acrocentric_p_spans(cytoband_path):
    """(chrom, 0, centromere_start) for the acrocentric short arms.

    These are satellite and rDNA repeat with no unique exonic sequence.
    The handful of bins that map there are repeat artefacts and they
    dominate the y-axis: on this cohort chr21p sits at -3.18 in the
    trisomy case and -2.91 in its mother, an identical excursion in both,
    which is how one can tell it is technical rather than biological.
    Conventional cytogenetics does not score these arms and neither does
    the caller. Leaving them in the figure would show a dramatic deletion
    that is not there, so they are excluded and the exclusion is stated
    on the figure rather than applied silently.
    """
    cb = load_cytoband(cytoband_path)
    spans = {}
    for chrom in ACROCENTRIC:
        cen = centromere(chrom, cb)
        if cen:
            spans[chrom] = (0, cen[0])
    return spans


def load(path, exclude_sex, drop_spans=None):
    """Point cloud per chromosome, from CNVkit .cnr or a GATK
    denoised-copy-ratio table.

    CNVkit columns: chromosome/start/end/depth/log2, and bins at zero
    depth are dropped -- their log2 is an artefact of dividing by the
    reference, not an observation.

    GATK columns: CONTIG/START/END/LOG2_COPY_RATIO, preceded by SAM-style
    '@' header lines. GATK has already dropped uninformative intervals at
    the FilterIntervals step, so there is no depth column to filter on.
    """
    chroms = [c for c in CHR_ORDER if not (exclude_sex and c in ("chrX", "chrY"))]
    pts = {c: [] for c in chroms}
    with open(path) as f:
        lines = [l for l in f if not l.startswith("@")]
    if not lines:
        return chroms, pts
    header = lines[0].rstrip("\n").split("\t")
    gatk = "LOG2_COPY_RATIO" in header
    cc, cs, ce, cl = (("CONTIG", "START", "END", "LOG2_COPY_RATIO") if gatk
                      else ("chromosome", "start", "end", "log2"))
    for r in csv.DictReader(lines, delimiter="\t"):
        c = r.get(cc)
        if c not in pts:
            continue
        try:
            l2 = float(r[cl])
        except (TypeError, ValueError, KeyError):
            continue
        if not gatk:
            try:
                if float(r["depth"]) <= 0:
                    continue
            except (TypeError, ValueError, KeyError):
                pass
        mid = (int(r[cs]) + int(r[ce])) / 2
        span = (drop_spans or {}).get(c)
        if span and span[0] <= mid < span[1]:
            continue
        pts[c].append((mid, l2))
    return chroms, pts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ratios", required=True,
                    help="CNVkit .cnr or GATK denoised copy-ratio table")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-chromosome-dir", default=None,
                    help="Also write one high-resolution panel per chromosome here. "
                         "The genome-wide trace compresses a chromosome into a few dozen "
                         "pixels, which is enough to see a whole-chromosome shift and not "
                         "enough to see one arm or a sub-arm segment. The per-chromosome "
                         "panels use a finer window for the same reason.")
    ap.add_argument("--per-chromosome-window-mb", type=float, default=0.25,
                    help="Window for the per-chromosome panels. Smaller than the "
                         "genome-wide window because the axis is ~50x wider per base, so "
                         "the trace can resolve a shorter event before noise dominates.")
    ap.add_argument("--window-mb", type=float, default=1.0)
    ap.add_argument("--min-bins-window", type=int, default=10)
    ap.add_argument("--ylim", type=float, default=2.0)
    ap.add_argument("--include-sex", action="store_true")
    ap.add_argument("--title", default=None)
    ap.add_argument("--cytoband", default=None,
                    help="Enables exclusion of the acrocentric short arms, which carry "
                         "no unique exonic sequence; see acrocentric_p_spans().")
    args = ap.parse_args()

    drop = acrocentric_p_spans(args.cytoband) if args.cytoband else {}
    chroms, pts = load(args.ratios, not args.include_sex, drop)
    offset, x0, ticks = {}, 0.0, []
    for c in chroms:
        offset[c] = x0
        ticks.append((x0 + CHR_LENGTHS[c] / 2e6, c.replace("chr", "")))
        x0 += CHR_LENGTHS[c] / 1e6
    total = x0

    # Baseline from window medians, not from raw bins: a handful of
    # saturated bins must not define where zero sits.
    win = {}
    for c in chroms:
        for mid, l2 in pts[c]:
            win.setdefault((c, int(mid // (args.window_mb * 1e6))), []).append(l2)
    wmed = {k: statistics.median(v) for k, v in win.items()
            if len(v) >= args.min_bins_window}
    baseline = statistics.median(wmed.values()) if wmed else 0.0
    noise = (1.4826 * statistics.median([abs(v - baseline) for v in wmed.values()])
             if wmed else float("nan"))

    fig, ax = plt.subplots(figsize=(22, 4.2))
    for c in chroms:
        xs = [offset[c] + m / 1e6 for m, _ in pts[c]]
        ys = [l2 - baseline for _, l2 in pts[c]]
        ax.scatter(xs, ys, s=0.7, c="#333333", alpha=0.18, linewidths=0, rasterized=True)
    for c in chroms:
        ks = sorted(k for k in wmed if k[0] == c)
        if not ks:
            continue
        xs = [offset[c] + (k[1] + 0.5) * args.window_mb for k in ks]
        ys = [wmed[k] - baseline for k in ks]
        ax.plot(xs, ys, color="#c0392b", linewidth=0.9, solid_joinstyle="round")

    for c in chroms:
        ax.axvline(offset[c], color="#bbbbbb", linewidth=0.6, zorder=0)
    ax.axhline(0, color="black", linewidth=0.7)
    for lvl, lab in ((0.585, "3/2"), (-1.0, "1/2")):
        ax.axhline(lvl, color="#7f8c8d", linewidth=0.5, linestyle=":")
        ax.text(total * 1.002, lvl, lab, fontsize=7, va="center", color="#7f8c8d")

    ax.set_xlim(0, total)
    ax.set_ylim(-args.ylim, args.ylim)
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], fontsize=8)
    ax.set_ylabel("Log$_2$ ratio", fontsize=10)
    ax.set_title(f"{args.title or args.sample}    "
                 f"[window {args.window_mb:g} Mb, noise MAD {noise:.3f}, "
                 f"single-copy gain = {0.585 / noise:.1f}x noise"
                 + ("; acrocentric short arms excluded]" if drop else "]"),
                 fontsize=10, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=160)
    print(f"{args.sample}: noise={noise:.4f} windows={len(wmed)} -> {args.out}")

    if args.per_chromosome_dir:
        os.makedirs(args.per_chromosome_dir, exist_ok=True)
        cb = load_cytoband(args.cytoband) if args.cytoband else {}
        pw = args.per_chromosome_window_mb
        for c in chroms:
            if not pts[c]:
                continue
            cw = {}
            for mid, l2 in pts[c]:
                cw.setdefault(int(mid // (pw * 1e6)), []).append(l2)
            cm = {k: statistics.median(v) for k, v in cw.items()
                  if len(v) >= max(4, args.min_bins_window // 3)}
            f2, a2 = plt.subplots(figsize=(16, 4.0))
            a2.scatter([m / 1e6 for m, _ in pts[c]], [l - baseline for _, l in pts[c]],
                       s=2.5, c="#333333", alpha=0.22, linewidths=0, rasterized=True)
            if cm:
                xs = sorted(cm)
                a2.plot([(k + 0.5) * pw for k in xs], [cm[k] - baseline for k in xs],
                        color="#c0392b", linewidth=1.3)
            # The centromere is drawn because an arm-level event is read
            # against it: a shift that stops at the centromere is a
            # different finding from one that crosses it.
            cen = centromere(c, cb) if cb else None
            if cen:
                a2.axvspan(cen[0] / 1e6, cen[1] / 1e6, color="#95a5a6", alpha=0.25, zorder=0)
            a2.axhline(0, color="black", linewidth=0.7)
            for lvl, lab in ((0.585, "3/2"), (-1.0, "1/2")):
                a2.axhline(lvl, color="#7f8c8d", linewidth=0.5, linestyle=":")
                a2.text(CHR_LENGTHS[c] / 1e6 * 1.003, lvl, lab, fontsize=8,
                        va="center", color="#7f8c8d")
            a2.set_xlim(0, CHR_LENGTHS[c] / 1e6)
            a2.set_ylim(-args.ylim, args.ylim)
            a2.set_xlabel("Position (Mb)", fontsize=10)
            a2.set_ylabel("Log$_2$ ratio", fontsize=10)
            a2.set_title(f"{args.title or args.sample} — {c}   "
                         f"[window {pw:g} Mb, genome noise MAD {noise:.3f}]",
                         fontsize=11, fontweight="bold")
            for sp in ("top", "right"):
                a2.spines[sp].set_visible(False)
            f2.tight_layout()
            f2.savefig(os.path.join(args.per_chromosome_dir, f"{c}.png"), dpi=170)
            plt.close(f2)
        print(f"{args.sample}: {len(chroms)} per-chromosome panels -> {args.per_chromosome_dir}")


if __name__ == "__main__":
    main()
