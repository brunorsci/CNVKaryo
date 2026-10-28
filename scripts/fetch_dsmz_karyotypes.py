#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collect published cytogenetic karyotypes for cell lines from the DSMZ
catalogue, via Cellosaurus for the accession lookup.

The karyotype is the orthogonal truth for this validation: it comes from
chromosome banding down a microscope, not from sequencing, so it is
independent of everything the pipeline does. DSMZ publishes a full ISCN
string plus the modal chromosome number for its leukemia-lymphoma lines.

The modal number matters and is captured separately. Cell lines are
frequently polyploid -- HL-60 is "human flat-moded hypotetraploid ...
82-88 XX, -X, -X, -8, -8, ..." -- so a '-8' there means two copies
against a tetraploid background, a relative loss. A read-depth method
measures everything against the sample's own baseline, so the modal
number is what converts a relative call into an absolute copy number.

Usage:
    python3 fetch_dsmz_karyotypes.py --cell-lines lines.txt --out dsmz.tsv
"""
import argparse
import csv
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request

CELLOSAURUS = "https://api.cellosaurus.org/search/cell-line?q={}&format=json&rows=5"
DSMZ = "https://www.dsmz.de/collection/catalogue/details/culture/{}"
UA = {"User-Agent": "cnvkaryo-validation/1.0 (research; contact via repository)"}
# DSMZ writes the modal number in three shapes, most commonly
# "46(45-48) XY", then "82-88 XX", then a bare "46,XX".
MODAL_PAREN = re.compile(r'(\d{2,3})\s*\(\s*(\d{2,3})\s*[-~]\s*(\d{2,3})\s*\)\s*,?\s*X')
MODAL = re.compile(r'(\d{2,3})\s*[-~]\s*(\d{2,3})\s*,?\s*X')
MODAL1 = re.compile(r'(?:^|[\s-])(\d{2,3})\s*,?\s*X')
PLOIDY_WORD = re.compile(
    r'(hypo|hyper|near|flat-moded)?\s*(diploid|triploid|tetraploid|pentaploid)', re.I)


def get(url, tries=3, pause=1.5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            if i == tries - 1:
                return None
            time.sleep(pause * (i + 1))
    return None


def dsmz_accession(name):
    """DSMZ ACC number for a cell line name, via Cellosaurus."""
    raw = get(CELLOSAURUS.format(urllib.parse.quote(name)))
    if not raw:
        return None, None
    try:
        entries = json.loads(raw)["Cellosaurus"]["cell-line-list"]
    except Exception:
        return None, None
    for entry in entries:
        for ref in entry.get("xref-list", []):
            if ref.get("database") == "DSMZ" and str(ref.get("accession", "")).startswith("ACC"):
                primary = entry.get("name-list", [{}])[0].get("value", name)
                return ref["accession"], primary
    return None, None


def dsmz_cytogenetics(acc):
    """The Cytogenetics free-text block from the DSMZ catalogue page."""
    page = get(DSMZ.format(acc))
    if not page:
        return None
    text = html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)))
    m = re.search(r"Cytogenetics:\s*(.+?)(?=\s(?:Viruses|References|Depositor|"
                  r"Isolation|Culture|Products|Relation|Comments|Authentication):|$)", text)
    return m.group(1).strip() if m else None


def modal_number(cyto):
    """(modal_low, modal_high, ploidy_word). ISCN writes the modal count
    before the sex chromosomes: '82-88 XX' or '46,XX'."""
    if not cyto:
        return None, None, None, None
    ploidy = None
    pm = PLOIDY_WORD.search(cyto)
    if pm:
        ploidy = " ".join(x for x in pm.groups() if x).strip().lower()
    m = MODAL_PAREN.search(cyto)
    if m:
        # "46(45-48)": 46 is the modal count, the parenthesis its range.
        return int(m.group(2)), int(m.group(3)), ploidy, int(m.group(1))
    m = MODAL.search(cyto)
    if m:
        lo, hi = int(m.group(1)), int(m.group(2))
        return lo, hi, ploidy, round((lo + hi) / 2)
    m = MODAL1.search(cyto)
    if m:
        return int(m.group(1)), int(m.group(1)), ploidy, int(m.group(1))
    return None, None, ploidy, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cell-lines", required=True,
                    help="One cell line name per line, or a TSV whose --name-column holds them")
    ap.add_argument("--name-column", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pause", type=float, default=1.0,
                    help="Seconds between requests, to stay polite to the catalogue.")
    args = ap.parse_args()

    names = []
    with open(args.cell_lines) as f:
        if args.name_column:
            for r in csv.DictReader(f, delimiter="\t"):
                v = (r.get(args.name_column) or "").strip()
                if v:
                    names.append(v)
        else:
            names = [l.strip() for l in f if l.strip()]

    fields = ["cell_line", "cellosaurus_name", "dsmz_accession", "modal_number",
              "modal_low", "modal_high", "ploidy", "karyotype"]
    rows, missing = [], []
    for i, name in enumerate(names, 1):
        acc, primary = dsmz_accession(name)
        cyto = dsmz_cytogenetics(acc) if acc else None
        lo, hi, ploidy, modal = modal_number(cyto)
        if not cyto:
            missing.append(name)
        rows.append(dict(cell_line=name, cellosaurus_name=primary or "",
                         dsmz_accession=acc or "", modal_number=modal if modal else "",
                         modal_low=lo if lo else "", modal_high=hi if hi else "",
                         ploidy=ploidy or "", karyotype=cyto or ""))
        print(f"[{i}/{len(names)}] {name:<16} {acc or '-':<8} "
              f"{'OK' if cyto else 'sem cariotipo'}", flush=True)
        time.sleep(args.pause)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows) - len(missing)}/{len(rows)} with a karyotype -> {args.out}")
    if missing:
        print("sem cariotipo: " + ", ".join(missing))


if __name__ == "__main__":
    main()
