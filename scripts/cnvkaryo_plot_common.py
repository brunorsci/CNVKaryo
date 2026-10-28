#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared drawing/classification module for cnvkaryo. Used by
plot_cnv_genome_scatter.py (grid + per-chromosome detail) and
export_cnv_tables.py (cn-to-band classification, cytogenetic naming) so the
grid, the per-chromosome detail, and the tables can never silently diverge.

Color/classification source of truth is always the absolute cn in
{sample}.call.cns (already recalibrated by cnvkit.py call, already
sex-aware via -x), never the raw .cns log2.
"""
import bisect

import csv

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# hg38 (GRCh38) chromosome lengths.
CHR_LENGTHS = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
}

GAIN_COLOR = "#1f77b4"
LOSS_COLOR = "#d62728"
NEUTRAL_COLOR = "#888888"
POINT_COLOR = "#555555"
GENE_LINE_COLOR = "#d4a017"


def load_gene_coords(path):
    """gene -> (chrom, start, end), from a TSV gene\\tchrom\\tstart\\tend."""
    coords = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 4:
                continue
            g, c, s, e = parts
            coords[g] = (c, int(s), int(e))
    return coords


def load_cytoband(path):
    """chrom -> sorted list of (start, end, band_name, gieStain), from the
    UCSC cytoBand format (chrom, start, end, name, gieStain)."""
    bands_by_chrom = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            stain = parts[4] if len(parts) > 4 else ""
            bands_by_chrom.setdefault(chrom, []).append((start, end, name, stain))
    for chrom in bands_by_chrom:
        bands_by_chrom[chrom].sort()
    return bands_by_chrom


def cytoband_label(chrom, start, end, cytobands):
    """Cytogenetic name of the [start,end) interval on chrom, e.g. '7q22.1'
    if it falls in a single band, or '7p22.3-7q36.3' if it spans several
    (uses the band touched by start and the one touched by end). Returns
    None if chrom has no cytoband data (e.g. alternate contigs)."""
    bands = cytobands.get(chrom)
    if not bands:
        return None
    starts = [b[0] for b in bands]

    def band_at(pos):
        idx = bisect.bisect_right(starts, pos) - 1
        if idx < 0:
            idx = 0
        if idx >= len(bands):
            idx = len(bands) - 1
        return bands[idx]

    first = band_at(start)
    last = band_at(max(start, end - 1))
    arm_prefix = chrom.replace("chr", "")
    if first[2] == last[2]:
        return f"{arm_prefix}{first[2]}"
    return f"{arm_prefix}{first[2]}-{arm_prefix}{last[2]}"


def ploidy_expected(chrom, sex):
    """Expected (neutral) copy number for this chromosome given the
    inferred sex ('Male'/'Female', as returned by cnvkit.py sex). Assumes
    a standard (non-male) CNVkit reference: autosome=2, chrX=2 if female
    or 1 if male, chrY=1 if male or 0 if female."""
    if chrom not in ("chrX", "chrY"):
        return 2
    sex_norm = (sex or "").strip().lower()
    is_male = sex_norm.startswith("m")
    if chrom == "chrX":
        return 1 if is_male else 2
    return 1 if is_male else 0


def classify_cn(cn, expected):
    """'gain' / 'loss' / 'neutral' from the absolute cn (integer, from
    call.cns) against the expected ploidy for that chromosome/sex."""
    if cn is None:
        return "neutral"
    if cn > expected:
        return "gain"
    if cn < expected:
        return "loss"
    return "neutral"


def classify_color(classification):
    return {"gain": GAIN_COLOR, "loss": LOSS_COLOR, "neutral": NEUTRAL_COLOR}.get(
        classification, NEUTRAL_COLOR
    )


def load_arm_calls(path):
    """chrom -> list of (start, end, call, iscn) for the regions the
    karyotype caller actually called. Read as text: pandas and csv both
    turn "+21" into a number given the chance, and the karyotype label
    then stops being a karyotype label."""
    calls = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["call"] not in ("gain", "loss"):
                continue
            key = "chromosome" if "chromosome" in r else "chrom"
            med = r.get("median_log2", "")
            calls.setdefault(r[key], []).append(
                (int(r["start"]), int(r["end"]), r["call"], r["iscn"],
                 float(med) if med not in ("", None) else 0.0))
    return calls


def render_chromosome_panel(ax, bands_chrom, cnr_chrom, chrom, y_min, y_max,
                             genes_here=None, point_size=4, seg_linewidth=2.5,
                             show_points=True, arm_calls_chrom=None):
    """Draws one panel (grid cell or standalone figure) for a single
    chromosome from the already-merged BAND table (not raw segments) --
    each band becomes exactly one colored line. .cnr points are background
    texture only and carry no classification weight.

    bands_chrom: list of dicts with at least start, end, classification,
                 log2 (the band's probe-weighted mean log2)
    cnr_chrom: list of (mid_mb, log2) .cnr bins for this chromosome, or
               None/[] to skip the background points
    genes_here: optional list of (gene, start, end) to highlight
    arm_calls_chrom: optional list of (start, end, call, iscn) from
                 call_arm_level.py. Shaded behind everything else so the
                 panel shows the karyotype verdict and the segment
                 evidence at once: the bands are the evidence, the shaded
                 span is what the tool concluded, and they are allowed to
                 disagree visibly.
    """
    for start, end, call, iscn, med in (arm_calls_chrom or []):
        colour = GAIN_COLOR if call == "gain" else LOSS_COLOR
        ax.axvspan(start / 1e6, end / 1e6, color=colour, alpha=0.13, zorder=0)
        # One line per call, drawn at the log2 the caller measured over the
        # whole region. The CBS bands underneath stay grey: inside a called
        # trisomy they include segments at +2.3 and at -1.9, neither of
        # which the tool asserts, and drawing them in gain/loss colour
        # makes the figure claim findings it rejected.
        ax.plot([start / 1e6, end / 1e6], [med, med], color=colour,
                linewidth=seg_linewidth * 1.6, solid_capstyle="butt", zorder=6)
        ax.annotate(iscn, xy=((start + end) / 2e6, y_max), xytext=(0, -2),
                    textcoords="offset points", ha="center", va="top",
                    fontsize=11, fontweight="bold", zorder=7, color=colour)
    if show_points and cnr_chrom:
        xs = [p[0] for p in cnr_chrom]
        ys = [p[1] for p in cnr_chrom]
        ax.scatter(xs, ys, s=point_size, c=POINT_COLOR, alpha=0.25, linewidths=0)

    # Every CBS band is grey. Colouring them by their own cn
    # classification painted alterations the tool never called -- on chr4
    # of a constitutional trisomy 21 a band sits at +0.67 across most of
    # the chromosome and read as a gain to anyone looking at the figure.
    # The bands are the evidence; the coloured line above is the finding.
    for band in bands_chrom:
        ax.plot([band["start"] / 1e6, band["end"] / 1e6],
                 [band["log2"], band["log2"]],
                 color=NEUTRAL_COLOR, linewidth=seg_linewidth, solid_capstyle="butt",
                 zorder=3, alpha=0.6)

    for gene, gs, ge in (genes_here or []):
        ax.axvspan(gs / 1e6, ge / 1e6, color=GENE_LINE_COLOR, alpha=0.4, zorder=1)

    ax.axhline(0, color="black", linewidth=0.6, zorder=2)
    ax.set_ylim(y_min, y_max)
