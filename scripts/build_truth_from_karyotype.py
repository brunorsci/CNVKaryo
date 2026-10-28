#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Turn published cell-line karyotypes into a copy-number truth table
expressed relative to each line's own ploidy.

A read-depth method has no absolute scale. It measures every bin against
the sample's own baseline, which is that sample's modal ploidy. A
karyotype must therefore be read the same way before it can be used as
truth.

This matters because cell lines are frequently polyploid. HL-60 is
"human flat-moded hypotetraploid ... 82-88 XX, -X, -X, -8, -8, -16, -17,
-17, +18, +22": the modal number near 85 puts its baseline at 4n, so
'-8, -8' means two copies of chromosome 8 against a background of four,
log2(2/4) = -1.0, a large loss. Read as though the line were diploid it
would be scored as a homozygous deletion, which is wrong, and in a
near-triploid line a single '+' is log2(4/3) = +0.42 rather than the
+0.585 of a diploid trisomy. Scoring every line against 2n would
misclassify roughly half of this cohort.

Base ploidy is taken as round(modal_number / 23) and the expected
relative log2 of chromosome N is log2(copies / base_ploidy), where
copies = base_ploidy + (gains of N) - (losses of N).

A cell line is a clone, so unlike a patient biopsy its published
karyotype describes the stemline rather than a composite of subclones --
there is no metaphase-fraction discount to apply. Sidelines, which DSMZ
reports after 'sideline with:', are a real subclonal fraction and are
excluded here rather than being counted as stemline events.

Usage:
    python3 build_truth_from_karyotype.py --karyotypes dsmz.tsv \
        --mapping sample_mapping.tsv --out-dir truth/
"""
import argparse
import csv
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnvkaryo_plot_common import CHR_LENGTHS, load_cytoband
from call_arm_level import centromere

GAIN = re.compile(r'(?<![\w(])\+(\d{1,2})(?![\w)])')
LOSS = re.compile(r'(?<![\w(])-(\d{1,2})(?![\w)])')
DEL = re.compile(r'del\((\d{1,2})\)\(([pq])', re.I)
DUP = re.compile(r'dup\((\d{1,2})\)\(([pq])', re.I)
SIDELINE = re.compile(r'\bsideline\b', re.I)

# Notation that leaves a chromosome's copy number undetermined by banding
# alone. A karyotype is a description of chromosome morphology, not of
# dosage: '-16, der(16)t(16;17)(q22;q22)' records the loss of a normal 16
# and the gain of a derivative that still carries part of 16, so the net
# copy number along that chromosome is not stated. Counting the '-16' as
# a whole-chromosome loss invents a truth the karyotype does not contain.
# Such chromosomes are marked indeterminate and excluded from scoring
# rather than guessed: 80% of the karyotypes in this cohort contain
# der(), 54% add() and 24% i().
DERIV = re.compile(r'\bder\((\d{1,2}|X|Y)[;)]', re.I)
ADD = re.compile(r'\badd\((\d{1,2}|X|Y)\)', re.I)
ISO = re.compile(r'\bi\((\d{1,2}|X|Y)\)', re.I)
INS = re.compile(r'\bins\(([\dXY;]+)\)', re.I)
TRANS_PARTNERS = re.compile(r'\b(?:t|der|ins)\(([\dXY;]+)\)', re.I)
MARKER = re.compile(r'\+\s*\d*\s*mar\b', re.I)
UNCERTAIN = re.compile(r'\?')
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def stemline(karyotype):
    """The stemline description, with any sideline text removed. DSMZ
    appends 'sideline with: ...' for a subclone; those events are present
    in only part of the culture and are not a fair test for bulk depth."""
    m = SIDELINE.search(karyotype)
    return karyotype[:m.start()] if m else karyotype


def indeterminate_chromosomes(text):
    """Chromosomes whose copy number this karyotype does not determine.

    Any chromosome appearing inside a derivative, an addition, an
    isochromosome or an insertion is included, together with every
    partner named in those brackets, because the exchanged material
    belongs to both. A reciprocal t(a;b) on its own is balanced and
    changes no dosage, so it does not make a chromosome indeterminate --
    it is in fact a useful negative: the karyotype says those
    chromosomes should read neutral by depth.
    """
    out = set()
    for rx in (DERIV, ADD, ISO):
        for m in rx.finditer(text):
            out.add(f"chr{m.group(1)}")
    for rx in (DERIV, INS):
        for m in rx.finditer(text):
            for part in re.split(r"[;)]", m.group(1)):
                part = part.strip()
                if part:
                    out.add(f"chr{part}")
    return {c for c in out if c in CHR_LENGTHS}


def chromosome_copies(text, base):
    """chrN -> copy number implied by the karyotype, relative to base."""
    copies = {c: base for c in AUTOSOMES}
    for n in GAIN.findall(text):
        c = f"chr{n}"
        if c in copies:
            copies[c] += 1
    for n in LOSS.findall(text):
        c = f"chr{n}"
        if c in copies:
            copies[c] -= 1
    return copies


def arm_events(text):
    """(chrom, arm, call) for explicit del()/dup() of a whole arm."""
    out = []
    for n, arm in DEL.findall(text):
        out.append((f"chr{n}", arm.lower(), "loss"))
    for n, arm in DUP.findall(text):
        out.append((f"chr{n}", arm.lower(), "gain"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--karyotypes", required=True)
    ap.add_argument("--mapping", required=True,
                    help="TSV with DNA_sample and Cell_line columns")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--cytoband", required=True,
                    help="Needed for the centromere, so an arm event is written over that "
                         "arm's real coordinates. Writing del(6)(q21) across the whole "
                         "chromosome makes 6p look like a documented loss as well.")
    ap.add_argument("--min-abs-log2", type=float, default=0.20,
                    help="An event whose expected relative log2 is smaller than this is "
                         "recorded but flagged not_resolvable: in a hyperpentaploid line "
                         "one extra copy of a chromosome is log2(6/5)=0.26, and in a "
                         "near-octaploid line log2(9/8)=0.17, below what read depth "
                         "separates from noise. Ploidy, not the tool, sets that limit.")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    cytobands = load_cytoband(args.cytoband)

    def arm_span(chrom, arm):
        cen = centromere(chrom, cytobands)
        if not cen:
            return None
        return (1, cen[0]) if arm == "p" else (cen[1], CHR_LENGTHS[chrom])
    line_to_sample = {}
    with open(args.mapping) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("Cell_line") and r.get("DNA_sample"):
                line_to_sample[r["Cell_line"]] = r["DNA_sample"]

    summary = []
    for r in csv.DictReader(open(args.karyotypes), delimiter="\t"):
        line = r["cell_line"]
        sample = line_to_sample.get(line)
        kar = r.get("karyotype", "")
        modal = r.get("modal_number", "")
        if not (sample and kar and modal):
            summary.append(dict(cell_line=line, sample=sample or "", base_ploidy="",
                                n_events="", status="no_karyotype_or_modal"))
            continue
        base = max(1, round(int(modal) / 23))
        stem = stemline(kar)
        copies = chromosome_copies(stem, base)
        unknown = indeterminate_chromosomes(stem)
        if MARKER.search(stem):
            # Marker chromosomes carry material of unidentified origin, so
            # no chromosome's total is pinned down by the karyotype.
            unknown |= set(AUTOSOMES)
        rows = []
        for chrom in sorted(unknown, key=lambda c: AUTOSOMES.index(c) if c in AUTOSOMES else 99):
            if chrom in AUTOSOMES:
                rows.append(dict(chrom=chrom, start=1, end=CHR_LENGTHS[chrom],
                                 expected_finding="indeterminate", region="whole",
                                 copies="", base_ploidy=base, expected_log2="",
                                 resolvable="n/a",
                                 source=f"DSMZ {r.get('dsmz_accession','')}: copy number "
                                        f"not determined by banding (derivative, addition, "
                                        f"isochromosome or marker involves this chromosome)"))
        for chrom in AUTOSOMES:
            cn = copies[chrom]
            if chrom in unknown or cn == base:
                continue
            exp = math.log2(cn / base) if cn > 0 else -5.0
            rows.append(dict(chrom=chrom, start=1, end=CHR_LENGTHS[chrom],
                             expected_finding="gain" if cn > base else "loss",
                             region="whole", copies=cn, base_ploidy=base,
                             expected_log2=round(exp, 3),
                             resolvable="yes" if abs(exp) >= args.min_abs_log2 else "no",
                             source=f"DSMZ {r.get('dsmz_accession','')} karyotype: {kar[:150]}"))
        for chrom, arm, call in arm_events(stem):
            if chrom not in CHR_LENGTHS or chrom in unknown:
                continue
            span = arm_span(chrom, arm)
            if not span:
                continue
            # A whole-chromosome count change on the same chromosome as a
            # partial deletion or duplication leaves the net dosage of that
            # arm unstated: EOL-1 is "+4 ... del(4)(q12)x2" and MOLM-13 is
            # "+8, +8, +8, del(8)(p1?p2?)". The karyotype records both the
            # extra copies and the partial loss without saying what the arm
            # ends up at, so it is marked indeterminate rather than scored
            # in one direction.
            finding = "indeterminate" if copies[chrom] != base else call
            rows.append(dict(chrom=chrom, start=span[0], end=span[1],
                             expected_finding=finding, region=arm, copies="",
                             base_ploidy=base, expected_log2="",
                             resolvable="n/a" if finding == "indeterminate" else "yes",
                             source=f"DSMZ {r.get('dsmz_accession','')} karyotype: {kar[:150]}"))
        fields = ["chrom", "start", "end", "expected_finding", "region", "copies",
                  "base_ploidy", "expected_log2", "resolvable", "source"]
        with open(os.path.join(args.out_dir, f"{sample}_truth.tsv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        summary.append(dict(cell_line=line, sample=sample, base_ploidy=base,
                            n_events=len(rows),
                            status="ok" if rows else "no_whole_or_arm_event"))

    with open(args.out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell_line", "sample", "base_ploidy",
                                          "n_events", "status"], delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    ok = [s for s in summary if s["status"] == "ok"]
    print(f"{len(ok)} lines with at least one scoreable event; "
          f"{len(summary) - len(ok)} without -> {args.out_summary}")


if __name__ == "__main__":
    main()
