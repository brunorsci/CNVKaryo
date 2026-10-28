#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Virtual karyotype: a 24-chromosome grid for one sample plus 24 individual
high-resolution per-chromosome images, in one pass (single read of the
.cnr/bands files).

Color/classification source: {sample}_cnv_bands.tsv (export_cnv_tables.py)
-- each band becomes exactly one colored line in the plot (not a raw .cns
segment). Individual .cnr points are background texture only and carry no
classification weight.

Usage:
    python3 plot_cnv_genome_scatter.py \
        --bands sample_cnv_bands.tsv --cnr sample.cnr \
        --output sample_genome_wide.png \
        --output-chromosomes-dir chromosomes/ \
        [--title "title"] [--gene-coords genes.tsv --gene GENE1 ...]

--gene/--gene-coords are optional and highlight genes from an
already-provided list (never "any gene" by default, only the ones passed
explicitly).
"""
import argparse
import csv
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from cnvkaryo_plot_common import (CHR_ORDER, CHR_LENGTHS, load_gene_coords,
                                  load_arm_calls, render_chromosome_panel)

Y_MIN_GRID, Y_MAX_GRID = -3, 3
Y_MIN_CHROM, Y_MAX_CHROM = -3, 3


def load_bands(path):
    bands_by_chrom = {c: [] for c in CHR_ORDER}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            chrom = r["chrom"]
            if chrom not in bands_by_chrom:
                continue
            bands_by_chrom[chrom].append({
                "start": int(r["start"]),
                "end": int(r["end"]),
                "classification": r["classification"],
                "log2": float(r["log2_weighted_mean"]),
            })
    return bands_by_chrom


def load_cnr_points(path):
    """chrom -> list of (mid_mb, log2), from a .cnr (chromosome,start,end,...,log2)."""
    points_by_chrom = {c: [] for c in CHR_ORDER}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            chrom = r["chromosome"]
            if chrom not in points_by_chrom:
                continue
            mid_mb = (int(r["start"]) + int(r["end"])) / 2 / 1e6
            points_by_chrom[chrom].append((mid_mb, float(r["log2"])))
    return points_by_chrom


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bands", required=True)
    ap.add_argument("--cnr", required=True)
    ap.add_argument("--output", required=True, help="24-chromosome grid image")
    ap.add_argument("--output-chromosomes-dir", required=True, help="Directory for the 24 individual images")
    ap.add_argument("--title", default=None)
    ap.add_argument("--gene-coords", default=None)
    ap.add_argument("--gene", action="append", default=[])
    ap.add_argument("--arms", default=None,
                    help="{sample}_arm_level.tsv. Shades the karyotype call behind the "
                         "bins and segments, so the panel shows the verdict and the "
                         "evidence it was read from at the same time.")
    args = ap.parse_args()

    os.makedirs(args.output_chromosomes_dir, exist_ok=True)

    bands_by_chrom = load_bands(args.bands)
    points_by_chrom = load_cnr_points(args.cnr)

    gene_coords = load_gene_coords(args.gene_coords) if args.gene_coords else {}
    genes_by_chrom = {}
    for g in args.gene:
        if g in gene_coords:
            c, s, e = gene_coords[g]
            genes_by_chrom.setdefault(c, []).append((g, s, e))
        else:
            sys.exit(f"Gene {g} not found in {args.gene_coords}")

    arm_calls = load_arm_calls(args.arms) if args.arms else {}

    # --- 4x6 grid ---
    ncols, nrows = 4, 6
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 22))
    axes_flat = axes.flatten()
    for idx, chrom in enumerate(CHR_ORDER):
        ax = axes_flat[idx]
        render_chromosome_panel(
            ax, bands_by_chrom[chrom], points_by_chrom[chrom], chrom,
            Y_MIN_GRID, Y_MAX_GRID, genes_here=genes_by_chrom.get(chrom, []),
            point_size=4, seg_linewidth=2.5,
            arm_calls_chrom=arm_calls.get(chrom, []),
        )
        ax.set_xlim(0, CHR_LENGTHS[chrom] / 1e6)
        ax.set_title(chrom.replace("chr", "Chr "), fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=7)
        if idx % ncols != 0:
            ax.set_yticklabels([])
    for idx in range(len(CHR_ORDER), nrows * ncols):
        axes_flat[idx].axis("off")

    legend_elements = [
        plt.Line2D([0], [0], color="#1f77b4", lw=3, label="Called gain"),
        plt.Line2D([0], [0], color="#d62728", lw=3, label="Called loss"),
        plt.Line2D([0], [0], color="#888888", lw=3, label="CBS segment, not called"),
        mpatches.Patch(color="#1f77b4", alpha=0.13,
                       label="shaded + thick line: the call, at its measured log2"),
    ]
    fig.legend(handles=legend_elements, loc="upper right", fontsize=11, bbox_to_anchor=(0.99, 0.995))
    fig.suptitle(args.title or "Virtual karyotype", fontsize=16, fontweight="bold", y=0.998)
    fig.supxlabel("Position (Mb)", fontsize=12)
    fig.supylabel("Copy ratio (log2)", fontsize=12)
    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.97])
    plt.savefig(args.output, dpi=130)
    plt.close(fig)
    print(f"Saved (grid): {args.output}")

    # --- 24 individual high-resolution images ---
    for chrom in CHR_ORDER:
        fig, ax = plt.subplots(figsize=(16, 5))
        render_chromosome_panel(
            ax, bands_by_chrom[chrom], points_by_chrom[chrom], chrom,
            Y_MIN_CHROM, Y_MAX_CHROM, genes_here=genes_by_chrom.get(chrom, []),
            point_size=6, seg_linewidth=3.5,
            arm_calls_chrom=arm_calls.get(chrom, []),
        )
        ax.set_xlim(0, CHR_LENGTHS[chrom] / 1e6)
        ax.set_xlabel("Position (Mb)")
        ax.set_ylabel("Copy ratio (log2)")
        ax.set_title(f"{args.title + ' -- ' if args.title else ''}{chrom}", fontsize=13, fontweight="bold")
        ax.legend(handles=legend_elements, loc="upper right", fontsize=9)
        plt.tight_layout()
        out_path = os.path.join(args.output_chromosomes_dir, f"{chrom}.png")
        plt.savefig(out_path, dpi=300)
        plt.close(fig)
    print(f"Saved (24 individual images): {args.output_chromosomes_dir}/")


if __name__ == "__main__":
    main()
