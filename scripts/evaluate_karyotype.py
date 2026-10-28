#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Score called karyotypes against documented truth, at chromosome-arm
resolution.

Every call is expanded to the arms it covers, so a whole-chromosome call
counts once per arm and is directly comparable to an arm-level one. A
documented region counts as an arm-level event when it covers at least
--min-coverage of that arm (the 50% convention CNApp attributes to
Davoli et al. 2017; published thresholds range from 5% to 80%, so it is
exposed as a parameter rather than hard-coded).

Acrocentric short arms are excluded: they carry no unique exonic
sequence and the caller never scores them. Sex chromosomes are excluded
unless --include-sex, matching the caller's default.

Outputs a per-sample breakdown and the pooled confusion matrix with
PPV, NPV, sensitivity and specificity.
"""
import argparse
import csv
import math
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnvkaryo_plot_common import CHR_ORDER, CHR_LENGTHS, load_cytoband
from call_arm_level import ACROCENTRIC, AUTOSOMES, centromere


def arm_intervals(cytobands, include_sex):
    chroms = AUTOSOMES + (["chrX", "chrY"] if include_sex else [])
    out = []
    for chrom in chroms:
        cen = centromere(chrom, cytobands)
        end = CHR_LENGTHS[chrom]
        if cen:
            if chrom not in ACROCENTRIC:
                out.append((chrom, "p", 0, cen[0]))
            out.append((chrom, "q", cen[1], end))
        else:
            out.append((chrom, "q", 0, end))
    return out


def truth_by_arm(truth_path, arms, min_coverage):
    """arm -> gain / loss / neutral / indeterminate.

    'indeterminate' marks a chromosome whose copy number the reference
    karyotype does not determine -- a derivative, addition,
    isochromosome or marker leaves the dosage unstated. Those arms are
    excluded from the confusion matrix entirely. Scoring them as neutral
    would manufacture false positives wherever the tool correctly
    detected unbalanced material, and scoring them as altered would
    manufacture false negatives; neither is a measurement.
    """
    t = {(c, a): "neutral" for c, a, _, _ in arms}
    if not truth_path or not os.path.exists(truth_path):
        return t
    regions = []
    with open(truth_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["expected_finding"] == "indeterminate":
                region = (r.get("region") or "whole").strip()
                for key in t:
                    if key[0] == r["chrom"] and (region == "whole" or key[1] == region):
                        t[key] = "indeterminate"
                continue
            if r["expected_finding"] in ("gain", "loss"):
                regions.append((r["chrom"], int(r["start"]), int(r["end"]),
                                r["expected_finding"]))
    for chrom, arm, s, e in arms:
        span = e - s
        if not span or t[(chrom, arm)] == "indeterminate":
            continue
        for label in ("gain", "loss"):
            covered = sum(max(0, min(e, re_) - max(s, rs))
                          for rc, rs, re_, rl in regions
                          if rc == chrom and rl == label)
            if covered / span >= min_coverage:
                t[(chrom, arm)] = label
    return t


def calls_by_arm(arms_path, arms):
    """Arm-level view of the caller's output: a whole-chromosome call is
    expanded onto every arm of that chromosome."""
    called = {(c, a): "neutral" for c, a, _, _ in arms}
    whole, per_arm = {}, {}
    with open(arms_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["call"] not in ("gain", "loss"):
                continue
            if r["region"] == "whole":
                whole[r["chrom"]] = r["call"]
            else:
                per_arm[(r["chrom"], r["region"])] = r["call"]
    for key in called:
        chrom, _ = key
        if chrom in whole:
            called[key] = whole[chrom]
        elif key in per_arm:
            called[key] = per_arm[key]
    return called


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--karyogram-dir", required=True,
                    help="Directory holding <sample>/<sample>_arm_level.tsv")
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--cytoband", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-coverage", type=float, default=0.5)
    ap.add_argument("--include-sex", action="store_true")
    ap.add_argument("--strata", default=None,
                    help="TSV (sample, stratum) splitting the panel by the kind of "
                         "karyotype being tested. Performance is not one number: a "
                         "constitutional aneuploidy and a subclonal composite cancer "
                         "karyotype are different problems for a depth-only method, "
                         "and pooling them hides both.")
    args = ap.parse_args()

    cytobands = load_cytoband(args.cytoband)
    arms = arm_intervals(cytobands, args.include_sex)

    samples = sorted(d for d in os.listdir(args.karyogram_dir)
                     if os.path.isfile(os.path.join(args.karyogram_dir, d,
                                                    f"{d}_arm_level.tsv")))
    strata = {}
    if args.strata:
        with open(args.strata) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                strata[r["sample"]] = r["stratum"]
    totals = dict(tp=0, fp=0, fn=0, tn=0)
    by_stratum = {}
    lines = []
    header = (f"{'sample':<40} {'TP':>3} {'FP':>3} {'FN':>3} {'TN':>4}  discrepancies")
    lines.append(header)
    lines.append("-" * len(header))
    for sample in samples:
        called = calls_by_arm(
            os.path.join(args.karyogram_dir, sample, f"{sample}_arm_level.tsv"), arms)
        truth_path = None
        for suffix in ("_truth.tsv", "_truth_array.tsv"):
            cand = os.path.join(args.truth_dir, f"{sample}{suffix}")
            if os.path.exists(cand):
                truth_path = cand
                break
        truth = truth_by_arm(truth_path, arms, args.min_coverage)
        s = dict(tp=0, fp=0, fn=0, tn=0, skipped=0)
        notes = []
        for key in sorted(called, key=lambda k: (CHR_ORDER.index(k[0]), k[1])):
            exp, got = truth[key], called[key]
            if exp == "indeterminate":
                s["skipped"] = s.get("skipped", 0) + 1
                continue
            tag = f"{key[0]}{key[1]}"
            if exp == "neutral" and got == "neutral":
                s["tn"] += 1
            elif exp == got:
                s["tp"] += 1
            elif exp == "neutral":
                s["fp"] += 1
                notes.append(f"FP {tag}={got}")
            elif got == "neutral":
                s["fn"] += 1
                notes.append(f"FN {tag} exp={exp}")
            else:
                s["fp"] += 1
                s["fn"] += 1
                notes.append(f"WRONG {tag} exp={exp} got={got}")
        for k in totals:
            totals[k] += s[k]
        st = strata.get(sample, "all")
        acc = by_stratum.setdefault(st, dict(tp=0, fp=0, fn=0, tn=0, n=0, skipped=0))
        acc["skipped"] += s.get("skipped", 0)
        for k in ("tp", "fp", "fn", "tn"):
            acc[k] += s[k]
        acc["n"] += 1
        src = os.path.basename(truth_path) if truth_path else "NO TRUTH FILE"
        lines.append(f"{sample:<40} {s['tp']:>3} {s['fp']:>3} {s['fn']:>3} {s['tn']:>4}  "
                     + (", ".join(notes) if notes else "-"))
        lines.append(f"{'':<40} truth: {src}")

    def wilson(k, n, z=1.96):
        """Wilson score interval. Reported because a proportion from a handful
        of events is not the same evidence as the same proportion from
        hundreds, and a point estimate alone hides that."""
        if not n:
            return None
        p = k / n
        d = 1 + z * z / n
        centre = (p + z * z / (2 * n)) / d
        half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
        return max(0.0, centre - half), min(1.0, centre + half)

    def block(title, c, n_samples):
        tp, fp, fn, tn = c["tp"], c["fp"], c["fn"], c["tn"]
        scored = tp + fp + fn + tn
        out = [f"{title} -- {n_samples} sample(s), {scored} arm decisions scored, "
               f"{c.get('skipped', 0)} skipped as indeterminate by the reference karyotype",
               f"  TP={tp}  FP={fp}  FN={fn}  TN={tn}"]
        for label, k, n in (("PPV (positive predictive value)", tp, tp + fp),
                            ("NPV (negative predictive value)", tn, tn + fn),
                            ("Sensitivity                    ", tp, tp + fn),
                            ("Specificity                    ", tn, tn + fp)):
            if not n:
                out.append(f"  {label} = n/a")
                continue
            lo, hi = wilson(k, n)
            out.append(f"  {label} = {k / n:.4f}   ({k}/{n}, 95% CI {lo:.3f}-{hi:.3f})")
        return out

    lines.append("")
    if by_stratum and set(by_stratum) != {"all"}:
        for st in sorted(by_stratum):
            c = by_stratum[st]
            lines.extend(block(f"STRATUM: {st}", c, c["n"]))
            lines.append("")
    lines.extend(block("POOLED (all strata)", totals, len(samples)))

    text = "\n".join(lines) + "\n"
    with open(args.out, "w") as f:
        f.write(text)
    print(text)


if __name__ == "__main__":
    main()
