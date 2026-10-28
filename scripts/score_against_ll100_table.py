#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Score cnvkaryo against the LL-100 publication's own sample table.

The table (Quentmeier et al. 2019, Sci Rep 9:8218) records one
Key_Alteration per line: the diagnostic lesion the authors attribute to
that cell line, for the material they sequenced. It is the right
reference for these runs, and it is not the same thing as the DSMZ
catalogue karyotype, which describes the collection's stock and can
differ by passage. Scoring against the catalogue produced two apparent
whole-chromosome failures that the publication does not claim.

Each key alteration is classed by what it does to DNA dosage:

  dosage      a deletion or amplification; read depth can detect it, so
              it counts as a positive and a miss is a false negative;
  balanced    a reciprocal translocation or fusion; it changes no copy
              number, so no depth method can see it. These are scored as
              negatives at the involved loci: calling dosage there would
              be a false positive. They are not scored as missed
              positives, because that would penalise the tool for a
              physical limit rather than a failure.
  amplicon    a fusion whose partner locus is also amplified in that
              line. Detectable, and scored as a positive, but only where
              the publication names the amplified partner.

Reporting positives and negatives separately, rather than pooling, is
what ACMG requires: a separate calculation per variant type.
"""
import argparse
import csv
import json
import os
import re
import sys

BALANCED = re.compile(r'^t\(|^inv\(|fusion|^[A-Z0-9]+-[A-Z0-9]+( e\d+-a\d+)?$', re.I)
DOSAGE = re.compile(r'\bdel\b|deletion|amp\b|amplif|gain|loss|biallelic|homozyg', re.I)


def classify(key_alteration):
    ka = (key_alteration or "").strip()
    if not ka or ka == "-":
        return "unannotated", []
    parts = [p.strip() for p in ka.split(",")]
    kinds = set()
    for p in parts:
        kinds.add("dosage" if DOSAGE.search(p) else "balanced")
    return ("mixed" if len(kinds) > 1 else kinds.pop()), parts


def genes_called(path, kinds=("homozygous_deletion", "focal_amplification")):
    """gene -> event_type, for unflagged focal calls."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("flag") or r.get("event_type") not in kinds:
                continue
            for g in (r.get("genes") or "").split(","):
                g = g.strip()
                if g and g != "Antitarget":
                    out[g] = r["event_type"]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-info", required=True)
    ap.add_argument("--karyogram-dir", required=True)
    ap.add_argument("--expected-genes", required=True,
                    help="JSON mapping a key alteration to the gene(s) that must be "
                         "called for it to count as recovered. Written by hand and "
                         "committed, so the mapping is auditable rather than inferred "
                         "from the free-text annotation at scoring time.")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    expected = json.load(open(args.expected_genes))
    info = {r["DNA_SRA"]: r for r in csv.DictReader(open(args.sample_info), delimiter="\t")}
    samples = sorted(d for d in os.listdir(args.karyogram_dir)
                     if os.path.isdir(os.path.join(args.karyogram_dir, d)))

    lines, tp, fn, scored_neg, fp_neg = [], 0, 0, 0, 0
    hdr = f"{'cell line':<13}{'class':<10}{'key alteration':<34}{'verdict'}"
    lines += [hdr, "-" * len(hdr)]
    for s in samples:
        r = info.get(s)
        if not r:
            continue
        cl = r["Cell_Line"]
        kind, parts = classify(r.get("Key_Alteration"))
        called = genes_called(os.path.join(args.karyogram_dir, s,
                                           f"{s}_homozygous_deletions.tsv"))
        verdicts = []
        for p in parts:
            want = expected.get(p)
            if want is None:
                verdicts.append(f"{p}: no gene mapping, not scored")
                continue
            hit = [g for g in want if g in called]
            if DOSAGE.search(p) or p in expected.get("_amplicon_keys", []):
                if hit:
                    tp += 1
                    verdicts.append(f"{p}: RECOVERED ({','.join(hit)})")
                else:
                    fn += 1
                    verdicts.append(f"{p}: MISSED (expected {'/'.join(want)})")
            else:
                scored_neg += 1
                if hit:
                    fp_neg += 1
                    verdicts.append(f"{p}: FALSE dosage call at a balanced locus")
                else:
                    verdicts.append(f"{p}: correctly silent (balanced)")
        if not parts:
            verdicts.append("no annotation, not scored")
        lines.append(f"{cl:<13}{kind:<10}{(r.get('Key_Alteration') or '-')[:33]:<34}"
                     f"{verdicts[0] if verdicts else ''}")
        for v in verdicts[1:]:
            lines.append(f"{'':<57}{v}")

    def ratio(a, b):
        return f"{a/b:.3f}" if b else "n/a"

    def wilson(k, n, z=1.96):
        if not n:
            return "n/a"
        import math
        p = k / n
        d = 1 + z * z / n
        c = (p + z * z / (2 * n)) / d
        h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
        return f"{max(0,c-h):.3f}-{min(1,c+h):.3f}"

    lines += ["", "Detectable (dosage) alterations named by the publication",
              f"  recovered {tp}, missed {fn}",
              f"  sensitivity = {ratio(tp, tp+fn)}   95% CI {wilson(tp, tp+fn)}",
              "",
              "Balanced alterations, scored as negatives at the involved loci",
              f"  loci scored {scored_neg}, false dosage calls {fp_neg}",
              f"  specificity = {ratio(scored_neg-fp_neg, scored_neg)}   "
              f"95% CI {wilson(scored_neg-fp_neg, scored_neg)}",
              "",
              "Balanced events are not counted as missed positives: a reciprocal",
              "translocation changes no dosage and no read-depth method can detect",
              "one. Counting them as failures would measure physics, not the tool."]
    text = "\n".join(lines) + "\n"
    open(args.out, "w").write(text)
    print(text)


if __name__ == "__main__":
    main()
