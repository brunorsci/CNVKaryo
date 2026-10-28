#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compares {sample}_cnv_bands.tsv against a published "truth" TSV (e.g. the
SEQC2/ascatNgs truth set for HCC1395, the expected whole-chr21 gain in a
trisomy 21, a published K562 CNV profile) and computes band-by-band
concordance. This is what backs a "validated" claim before any
submission -- not just a visual check.

Expected truth.tsv (hand-curated from the literature/consortium for each
validation sample): chrom, start, end, expected_finding (gain/loss/
neutral), source

Usage:
    python3 compare_bands_to_truth.py --bands sample_cnv_bands.tsv \
        --truth truth.tsv --sample SAMPLE --output sample_concordance.tsv
"""
import argparse
import csv


def load_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def overlap_bp(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bands", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    bands = load_tsv(args.bands)
    for b in bands:
        b["start"] = int(b["start"])
        b["end"] = int(b["end"])

    truth = load_tsv(args.truth)
    for t in truth:
        t["start"] = int(t["start"])
        t["end"] = int(t["end"])

    rows = []
    matched_truth = 0
    for t in truth:
        overlapping = [
            b for b in bands
            if b["chrom"] == t["chrom"] and overlap_bp(b["start"], b["end"], t["start"], t["end"]) > 0
        ]
        total_overlap = sum(overlap_bp(b["start"], b["end"], t["start"], t["end"]) for b in overlapping)
        truth_len = t["end"] - t["start"]
        overlap_frac = total_overlap / truth_len if truth_len > 0 else 0.0
        classifications = sorted({b["classification"] for b in overlapping})
        matched = t["expected_finding"] in classifications
        if matched:
            matched_truth += 1
        rows.append({
            "sample": args.sample,
            "chrom": t["chrom"],
            "truth_start": t["start"],
            "truth_end": t["end"],
            "expected_finding": t["expected_finding"],
            "source": t.get("source", ""),
            "overlapping_bands": ";".join(b["band_id"] for b in overlapping),
            "classifications_found": ";".join(classifications),
            "overlap_bp": total_overlap,
            "overlap_fraction_of_truth_region": round(overlap_frac, 4),
            "matched": matched,
        })

    truth_band_ids = set()
    for t in truth:
        for b in bands:
            if b["chrom"] == t["chrom"] and overlap_bp(b["start"], b["end"], t["start"], t["end"]) > 0:
                truth_band_ids.add(b["band_id"])
    unmatched_bands = [b for b in bands if b["band_id"] not in truth_band_ids and b["classification"] != "neutral"]

    header = ["sample", "chrom", "truth_start", "truth_end", "expected_finding", "source",
              "overlapping_bands", "classifications_found", "overlap_bp",
              "overlap_fraction_of_truth_region", "matched"]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    pct = 100 * matched_truth / len(truth) if truth else 0.0
    print(f"[{args.sample}] {matched_truth}/{len(truth)} expected findings correctly flagged ({pct:.1f}%)")
    print(f"[{args.sample}] {len(unmatched_bands)} of our altered band(s) have no match in the published truth "
          f"(may be a novel finding, or a region the truth set doesn't cover -- not necessarily an error):")
    for b in unmatched_bands:
        print(f"    {b['band_id']} {b['chrom']}:{b['start']}-{b['end']} ({b['classification']}, {b.get('cytoband', '')})")
    print(f"  -> {args.output}")


if __name__ == "__main__":
    main()
