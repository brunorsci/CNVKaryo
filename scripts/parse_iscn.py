#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse an ISCN karyotype string into per-event clonal fractions.

A documented karyotype is not by itself a fair test for a depth-based
method. ISCN records what was seen down a microscope, metaphase by
metaphase, and a composite karyotype pools clones that are present in
different cells. `46,XX,t(15;17)[18]/47,idem,+8[2]` documents a +8, but
in 2 of 20 metaphases -- roughly 10% of cells, which shifts bulk read
depth by about 0.07 log2 and cannot be detected. Scoring such an event
as a miss measures the sample, not the tool.

This converts the string into events with the fraction of scored
metaphases carrying each one, so a validation panel can be built from
clonal events and the subclonal ones reported separately rather than
silently counted as false negatives.

Clone syntax handled:
  - clones separated by '/', each ending in '[n]' scored metaphases
  - 'idem', 'sl', 'sdl' meaning "the previous clone, plus what follows"
  - '[cp n]' composite clones, where ISCN itself states the events were
    pooled across cells; these are marked composite and their fraction
    is an upper bound, not a measurement
  - ranges such as '45~50' or '44-46' in the modal number, which carry
    no clone information and are ignored

Usage:
    python3 parse_iscn.py --karyotype "47,XY,+21[20]"
    python3 parse_iscn.py --tsv samples.tsv --karyotype-column karyotype \
        --out events.tsv
"""
import argparse
import csv
import re
import sys

GAIN = re.compile(r'(?<![\w(])\+(\d{1,2})(?![\w)])')
LOSS = re.compile(r'(?<![\w(])-(\d{1,2})(?![\w)])')
SEX_LOSS = re.compile(r'(?<![\w(])-([XY])(?![\w)])')
DEL = re.compile(r'del\((\d{1,2})\)\(([pq])', re.I)
DUP = re.compile(r'dup\((\d{1,2})\)\(([pq])', re.I)
ADD = re.compile(r'add\((\d{1,2})\)\(([pq])', re.I)
COUNT = re.compile(r'\[\s*(cp)?\s*(\d+)\s*\]', re.I)
IDEM = re.compile(r'\b(idem|sl|sdl\d*)\b', re.I)


def split_clones(karyotype):
    """[(clone_text, n_metaphases, is_composite)] in document order."""
    clones = []
    for part in karyotype.split("/"):
        part = part.strip()
        if not part:
            continue
        counts = COUNT.findall(part)
        if counts:
            is_cp = any(c[0] for c in counts)
            n = max(int(c[1]) for c in counts)
        else:
            is_cp, n = False, None
        clones.append((COUNT.sub("", part).strip(), n, is_cp))
    return clones


def events_in(text):
    """Whole-chromosome and arm-level events written in this clone."""
    out = set()
    for c in GAIN.findall(text):
        out.add(("gain", f"chr{c}", "whole"))
    for c in LOSS.findall(text):
        out.add(("loss", f"chr{c}", "whole"))
    for c in SEX_LOSS.findall(text):
        out.add(("loss", f"chr{c.upper()}", "whole"))
    for c, arm in DEL.findall(text):
        out.add(("loss", f"chr{c}", arm.lower()))
    for c, arm in DUP.findall(text):
        out.add(("gain", f"chr{c}", arm.lower()))
    return out


def parse(karyotype):
    """[{call, chrom, region, metaphases, total_metaphases, fraction,
    composite}] for one ISCN string."""
    clones = split_clones(karyotype)
    if not clones:
        return [], 0
    # ISCN omits the scored count when only one clone was reported; treat
    # that as a single fully-represented clone rather than dropping it.
    if all(n is None for _, n, _ in clones):
        clones = [(t, 1, cp) for t, n, cp in clones]
    clones = [(t, (n if n is not None else 0), cp) for t, n, cp in clones]
    total = sum(n for _, n, _ in clones)
    if not total:
        return [], 0

    carried = {}
    inherited = set()
    for text, n, is_cp in clones:
        own = events_in(text)
        # 'idem'/'sl' means this clone also carries the previous clone's
        # events; without expanding it a subclone's inherited events are
        # invisible and their fraction comes out too low.
        evs = (inherited | own) if IDEM.search(text) else own
        inherited = evs
        for e in evs:
            rec = carried.setdefault(e, {"n": 0, "composite": False})
            rec["n"] += n
            rec["composite"] = rec["composite"] or is_cp

    out = []
    for (call, chrom, region), rec in sorted(carried.items()):
        out.append(dict(call=call, chrom=chrom, region=region,
                        metaphases=rec["n"], total_metaphases=total,
                        fraction=round(rec["n"] / total, 3),
                        composite=rec["composite"]))
    return out, total


def classify(fraction, composite, clonal_min, detectable_min):
    """How a validation panel should treat this event."""
    if composite:
        return "composite_upper_bound"
    if fraction >= clonal_min:
        return "clonal"
    if fraction >= detectable_min:
        return "subclonal_borderline"
    return "subclonal_undetectable"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--karyotype")
    ap.add_argument("--tsv")
    ap.add_argument("--karyotype-column", default="karyotype")
    ap.add_argument("--id-column", default="biosample")
    ap.add_argument("--out")
    ap.add_argument("--clonal-min", type=float, default=0.80,
                    help="At or above this fraction of metaphases an event is treated as "
                         "clonal and a fair test for a bulk depth method. A clonal "
                         "single-copy gain shifts log2 by ~0.585 at fraction 1.0 and "
                         "~0.49 at 0.80, both far above this pipeline's measured "
                         "negative-control floor of 0.131.")
    ap.add_argument("--detectable-min", type=float, default=0.40,
                    help="Below this the expected shift falls under the noise floor and "
                         "the event is not a fair test. At fraction 0.40 a single-copy "
                         "gain shifts log2 by ~0.26, already at the calling threshold.")
    args = ap.parse_args()

    fields = ["id", "karyotype", "call", "chrom", "region", "metaphases",
              "total_metaphases", "fraction", "composite", "panel_class"]
    rows = []
    if args.karyotype:
        evs, _ = parse(args.karyotype)
        for e in evs:
            rows.append(dict(id="-", karyotype=args.karyotype, **e,
                             panel_class=classify(e["fraction"], e["composite"],
                                                  args.clonal_min, args.detectable_min)))
    elif args.tsv:
        with open(args.tsv) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                k = (r.get(args.karyotype_column) or "").strip()
                if not k:
                    continue
                evs, _ = parse(k)
                for e in evs:
                    rows.append(dict(id=r.get(args.id_column, ""), karyotype=k, **e,
                                     panel_class=classify(e["fraction"], e["composite"],
                                                          args.clonal_min,
                                                          args.detectable_min)))
    else:
        sys.exit("Need --karyotype or --tsv")

    out = open(args.out, "w", newline="") if args.out else sys.stdout
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
    if args.out:
        out.close()
        print(f"{len(rows)} events -> {args.out}")


if __name__ == "__main__":
    main()
