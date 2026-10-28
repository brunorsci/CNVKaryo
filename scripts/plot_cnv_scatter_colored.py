#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Whole-chromosome CNV scatter, in the style of `cnvkit.py scatter`, with two
differences: (1) the band is colored by state -- blue = gain, red = loss,
gray = neutral -- from {sample}_cnv_bands.tsv (absolute cn, sex-aware, not
an arbitrary log2 cutoff here); (2) the gene label is anchored at the
height of the band that actually overlaps the gene.

Only plots genes passed explicitly via --gene (never "any gene" by
default) -- a gene must be present in --gene-coords, otherwise this aborts
with a clear error rather than silently producing a plot with no
highlight.

Usage:
    python3 plot_cnv_scatter_colored.py \
        --cnr sample.cnr --bands sample_cnv_bands.tsv --chrom chr18 \
        --gene BCL2 --gene-coords gene_coords_all.tsv \
        --output output.png [--title "custom title"]

Accepts more than one --gene (all must be on the same --chrom).
"""
import argparse
import csv
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cnvkaryo_plot_common import load_gene_coords, classify_color

Y_MIN, Y_MAX = -5, 5  # same automatic cutoff as cnvkit.py scatter (excludes extreme outliers like -20, which are low-count noise, not linear signal)


def load_bands_chrom(path, chrom):
    bands = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r["chrom"] != chrom:
                continue
            bands.append({
                "start": int(r["start"]),
                "end": int(r["end"]),
                "log2": float(r["log2_weighted_mean"]),
                "classification": r["classification"],
            })
    return bands


def load_cnr_chrom(path, chrom):
    points = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r["chromosome"] != chrom:
                continue
            mid = (int(r["start"]) + int(r["end"])) / 2 / 1e6
            points.append((mid, float(r["log2"])))
    return points


def label_y_for_gene(bands_chrom, gene_start, gene_end):
    """log2 of the band that overlaps the gene (the one with the largest
    overlap, if more than one) -- anchors the label at the right height."""
    best, best_overlap = None, -1
    for band in bands_chrom:
        overlap = min(band["end"], gene_end) - max(band["start"], gene_start)
        if overlap > best_overlap:
            best_overlap, best = overlap, band
    if best is None or best_overlap <= 0:
        return None
    return best["log2"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cnr", required=True, help=".cnr file (log2 per bin)")
    ap.add_argument("--bands", required=True, help="{sample}_cnv_bands.tsv (export_cnv_tables.py)")
    ap.add_argument("--chrom", required=True, help="Chromosome to plot (e.g. chr18)")
    ap.add_argument("--gene", action="append", required=True, help="Gene(s) of interest to highlight (repeat for more than one)")
    ap.add_argument("--gene-coords", required=True, help="TSV gene\\tchrom\\tstart\\tend")
    ap.add_argument("--output", required=True)
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    coords = load_gene_coords(args.gene_coords)
    for gene in args.gene:
        if gene not in coords:
            sys.exit(f"ERROR: gene {gene!r} not found in {args.gene_coords} -- "
                     f"check the symbol (it must match the coordinates table)")

    cnr_points = load_cnr_chrom(args.cnr, args.chrom)
    bands_chrom = load_bands_chrom(args.bands, args.chrom)

    fig, ax = plt.subplots(figsize=(10, 7))

    if cnr_points:
        xs = [p[0] for p in cnr_points]
        ys = [p[1] for p in cnr_points]
        ax.scatter(xs, ys, s=8, c="#555555", alpha=0.3, linewidths=0)

    for band in bands_chrom:
        ax.plot([band["start"] / 1e6, band["end"] / 1e6], [band["log2"], band["log2"]],
                color=classify_color(band["classification"]), linewidth=3,
                solid_capstyle="butt", zorder=5)

    plotted_genes = []
    for gene in args.gene:
        g_chrom, gs, ge = coords[gene]
        if g_chrom != args.chrom:
            sys.exit(f"ERROR: gene {gene!r} is on {g_chrom}, not {args.chrom} (--chrom)")
        ax.axvspan(gs / 1e6, ge / 1e6, color="#d4a017", alpha=0.4, zorder=1)

        anchor_y = label_y_for_gene(bands_chrom, gs, ge)
        if anchor_y is None:
            anchor_y = 0.0  # no overlapping band (rare) -- falls back to neutral
        label_y = anchor_y - 0.5 if anchor_y > Y_MAX - 1 else anchor_y + 0.3
        label_y = max(Y_MIN + 0.2, min(Y_MAX - 0.2, label_y))
        ax.text((gs + ge) / 2 / 1e6, label_y, gene, ha="center", fontsize=10,
                fontweight="bold", zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.75))
        plotted_genes.append(gene)

    ax.axhline(0, color="black", linewidth=0.8, zorder=2)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_xlabel("Position (Mb)")
    ax.set_ylabel("Copy ratio (log2)")
    ax.set_title(args.title or f"{args.chrom} -- {', '.join(plotted_genes)}")

    legend_elements = [
        plt.Line2D([0], [0], color="#1f77b4", lw=3, label="Gain (cn > expected ploidy)"),
        plt.Line2D([0], [0], color="#d62728", lw=3, label="Loss (cn < expected ploidy)"),
        plt.Line2D([0], [0], color="#888888", lw=3, label="Neutral"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
