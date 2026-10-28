#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Produces the 3 audit tables for the virtual karyotype, from
{sample}.call.cns (already recalibrated and sex-aware via cnvkit.py call):

  {sample}_cnv_segments.tsv     -- one row per raw call.cns segment
  {sample}_cnv_bands.tsv        -- adjacent segments with the SAME
                                   classification (gain/loss/neutral)
                                   collapsed into one band -- this is what
                                   the plots draw: one table row = one
                                   colored bar in the image
  {sample}_cnv_breakpoints.tsv  -- every band-to-band boundary where the
                                   classification changes. Does not
                                   distinguish a balanced rearrangement
                                   from a simple deletion/duplication edge
                                   -- it only flags where copy number
                                   changes.

Usage:
    python3 export_cnv_tables.py --call-cns sample.call.cns --sex-file sample_sex.txt
        --cytoband cytoBand_hg38.txt --sample SAMPLE
        --out-segments out_segments.tsv --out-bands out_bands.tsv
        --out-breakpoints out_breakpoints.tsv
"""
import argparse
import csv
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnvkaryo_plot_common import (
    CHR_ORDER, load_cytoband, cytoband_label, ploidy_expected, classify_cn,
)

SEGMENTS_HEADER = [
    "sample", "band_id", "chrom", "start", "end", "length_bp", "log2", "cn",
    "depth", "probes", "weight", "p_ttest", "genes",
]
BANDS_HEADER = [
    "sample", "band_id", "chrom", "start", "end", "length_bp", "cytoband",
    "classification", "n_segments", "cn_min", "cn_max",
    "log2_weighted_mean", "p_ttest_most_conservative", "probes_total",
    "genes_union",
]
BREAKPOINTS_HEADER = [
    "sample", "chrom", "boundary_start", "boundary_end", "cytoband",
    "classification_before", "classification_after", "cn_before", "cn_after",
    "log2_before", "log2_after", "band_id_before", "band_id_after",
]
BREAKPOINTS_NOTE = (
    "# NOTE: this table only flags where copy number changes (the boundary "
    "between bands of different classification). A BALANCED rearrangement "
    "(reciprocal translocation/inversion) carries no read-depth signature -- "
    "it cannot be distinguished from a plain deletion/duplication edge using "
    "this data alone. Findings here need orthogonal validation (e.g. "
    "conventional karyotype, FISH) for full characterization.\n"
)


def parse_sex(sex_file):
    with open(sex_file) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    if len(lines) < 2:
        sys.exit(f"Sex file {sex_file} is missing the expected data row")
    header = lines[0].split("\t")
    row = lines[-1].split("\t")
    rec = dict(zip(header, row))
    sex = rec.get("sex")
    if sex not in ("Male", "Female"):
        sys.exit(f"Invalid inferred sex in {sex_file}: {sex!r}")
    return sex


def load_call_cns(path):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            chrom = r["chromosome"]
            if chrom not in CHR_ORDER:
                continue  # drop alt/random/Un contigs -- no cytoband, no statistical power
            rows.append({
                "chrom": chrom,
                "start": int(r["start"]),
                "end": int(r["end"]),
                "gene": r.get("gene", ""),
                "log2": float(r["log2"]),
                "cn": int(r["cn"]),
                "depth": float(r.get("depth", "nan")) if r.get("depth") not in (None, "") else None,
                "p_ttest": float(r["p_ttest"]) if r.get("p_ttest") not in (None, "") else None,
                "probes": int(r["probes"]) if r.get("probes") not in (None, "") else None,
                "weight": float(r["weight"]) if r.get("weight") not in (None, "") else None,
            })
    rows.sort(key=lambda r: (CHR_ORDER.index(r["chrom"]), r["start"]))
    return rows


def build_tables(segments, sex, cytobands):
    segments_out = []
    bands_out = []
    breakpoints_out = []

    by_chrom = {}
    for seg in segments:
        by_chrom.setdefault(seg["chrom"], []).append(seg)

    for chrom in CHR_ORDER:
        chrom_segs = by_chrom.get(chrom, [])
        if not chrom_segs:
            continue
        expected = ploidy_expected(chrom, sex)
        for seg in chrom_segs:
            seg["classification"] = classify_cn(seg["cn"], expected)

        band_idx = 0
        current = None
        chrom_bands = []
        for seg in chrom_segs:
            if current is not None and seg["classification"] == current["classification"]:
                current["segments"].append(seg)
                current["end"] = seg["end"]
            else:
                if current is not None:
                    chrom_bands.append(current)
                band_idx += 1
                current = {
                    "band_id": f"{chrom}_b{band_idx}",
                    "chrom": chrom,
                    "start": seg["start"],
                    "end": seg["end"],
                    "classification": seg["classification"],
                    "segments": [seg],
                }
        if current is not None:
            chrom_bands.append(current)

        for band in chrom_bands:
            segs = band["segments"]
            probes_total = sum((s["probes"] or 0) for s in segs)
            if probes_total > 0:
                log2_w = sum(s["log2"] * (s["probes"] or 0) for s in segs) / probes_total
            else:
                log2_w = sum(s["log2"] for s in segs) / len(segs)
            p_values = [s["p_ttest"] for s in segs if s["p_ttest"] is not None]
            p_worst = max(p_values) if p_values else None
            genes_union = sorted({
                g.strip() for s in segs for g in (s["gene"] or "").split(",") if g.strip()
            })
            cyto = cytoband_label(chrom, band["start"], band["end"], cytobands)

            for seg in segs:
                segments_out.append({
                    "band_id": band["band_id"],
                    "chrom": chrom,
                    "start": seg["start"],
                    "end": seg["end"],
                    "length_bp": seg["end"] - seg["start"],
                    "log2": seg["log2"],
                    "cn": seg["cn"],
                    "depth": seg["depth"],
                    "probes": seg["probes"],
                    "weight": seg["weight"],
                    "p_ttest": seg["p_ttest"],
                    "genes": seg["gene"],
                })

            bands_out.append({
                "band_id": band["band_id"],
                "chrom": chrom,
                "start": band["start"],
                "end": band["end"],
                "length_bp": band["end"] - band["start"],
                "cytoband": cyto or "",
                "classification": band["classification"],
                "n_segments": len(segs),
                "cn_min": min(s["cn"] for s in segs),
                "cn_max": max(s["cn"] for s in segs),
                "log2_weighted_mean": log2_w,
                "p_ttest_most_conservative": p_worst,
                "probes_total": probes_total,
                "genes_union": ";".join(genes_union),
            })

        for i in range(len(chrom_bands) - 1):
            a, b = chrom_bands[i], chrom_bands[i + 1]
            if a["classification"] == b["classification"]:
                continue
            seg_a, seg_b = a["segments"][-1], b["segments"][0]
            breakpoints_out.append({
                "chrom": chrom,
                "boundary_start": a["end"],
                "boundary_end": b["start"],
                "cytoband": cytoband_label(chrom, a["end"], max(a["end"], b["start"]), cytobands) or "",
                "classification_before": a["classification"],
                "classification_after": b["classification"],
                "cn_before": seg_a["cn"],
                "cn_after": seg_b["cn"],
                "log2_before": seg_a["log2"],
                "log2_after": seg_b["log2"],
                "band_id_before": a["band_id"],
                "band_id_after": b["band_id"],
            })

    return segments_out, bands_out, breakpoints_out


def write_tsv(path, header, rows, sample, note=None):
    with open(path, "w", newline="") as f:
        if note:
            f.write(note)
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            row = dict(row)
            row["sample"] = sample
            writer.writerow(row)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--call-cns", required=True)
    ap.add_argument("--sex-file", required=True)
    ap.add_argument("--declared-sex", default=None,
                    help="Overrides --sex-file. Must match what cnvkit.py call was given, "
                         "otherwise chrX/chrY expected ploidy here disagrees with the cn it calls.")
    ap.add_argument("--cytoband", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-segments", required=True)
    ap.add_argument("--out-bands", required=True)
    ap.add_argument("--out-breakpoints", required=True)
    args = ap.parse_args()

    sex = args.declared_sex or parse_sex(args.sex_file)
    if sex not in ("Male", "Female"):
        sys.exit(f"Invalid sex {sex!r}: must be 'Male' or 'Female'")
    cytobands = load_cytoband(args.cytoband)
    segments = load_call_cns(args.call_cns)
    segments_out, bands_out, breakpoints_out = build_tables(segments, sex, cytobands)

    write_tsv(args.out_segments, SEGMENTS_HEADER, segments_out, args.sample)
    write_tsv(args.out_bands, BANDS_HEADER, bands_out, args.sample)
    write_tsv(args.out_breakpoints, BREAKPOINTS_HEADER, breakpoints_out, args.sample, note=BREAKPOINTS_NOTE)

    print(f"[{args.sample}] sex={sex} segments={len(segments_out)} bands={len(bands_out)} breakpoints={len(breakpoints_out)}")
    print(f"  -> {args.out_segments}")
    print(f"  -> {args.out_bands}")
    print(f"  -> {args.out_breakpoints}")


if __name__ == "__main__":
    main()
