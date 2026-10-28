#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove correlated technical bias from exome log2 ratios by projecting out
the leading principal components of a panel of samples.

A pooled reference corrects each bin independently: it subtracts that
bin's median across the panel. That removes bias which is constant per
bin, but not bias which is *shared across bins and varies between
samples* -- capture efficiency drifting with library prep, GC response
differing by batch, insert-size effects. Those make whole groups of bins
move together, and a per-bin median cannot see them.

Measured here: against a flat reference the median window noise over ten
LL-100 exomes was 0.1413 log2; a pooled reference over the same 100
samples brought it only to 0.1026, short of what a legible figure needs.
The residual is correlated structure, which is what this removes.

This is the method behind GATK's DenoiseReadCounts
(CreateReadCountPanelOfNormals performs an SVD of the panel and
DenoiseReadCounts projects the case sample onto it); it is implemented
directly here rather than adding a dependency, and the panel is built
from the same coverage files the rest of the pipeline already produces.

The case sample must not be part of the panel. Removing components
derived from a matrix containing the sample would also remove that
sample's real copy number, which is the one thing that must survive.

Usage:
    python3 denoise_svd.py --panel panel/*.cnr --case sample.cnr \
        --out sample.denoised.cnr --components 10
"""
import argparse
import csv
import os
import sys

import numpy as np


def read_cnr(path):
    keys, log2, rows = [], [], []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            keys.append((r["chromosome"], r["start"], r["end"]))
            log2.append(float(r["log2"]))
            rows.append(r)
    return keys, np.asarray(log2, dtype=float), rows


def build_panel(paths, keys):
    """bins x samples matrix, aligned to the case sample's bins."""
    index = {k: i for i, k in enumerate(keys)}
    cols = []
    for p in paths:
        k, v, _ = read_cnr(p)
        col = np.full(len(keys), np.nan)
        for kk, vv in zip(k, v):
            i = index.get(kk)
            if i is not None:
                col[i] = vv
        cols.append(col)
    return np.column_stack(cols)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", nargs="+", required=True)
    ap.add_argument("--case", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--components", type=int, default=10,
                    help="Number of leading components to project out. Too few leaves "
                         "bias; too many starts removing real copy number, since a "
                         "component can align with an event shared by several panel "
                         "members. Chosen by measuring window noise, not by default.")
    ap.add_argument("--clip", type=float, default=5.0,
                    help="Clip |log2| before the decomposition so a handful of saturated "
                         "bins cannot define a component.")
    args = ap.parse_args()

    panel_paths = [p for p in args.panel if os.path.abspath(p) != os.path.abspath(args.case)]
    if len(panel_paths) < args.components + 2:
        sys.exit(f"Panel of {len(panel_paths)} is too small for {args.components} components")

    keys, case, rows = read_cnr(args.case)
    M = build_panel(panel_paths, keys)

    usable = ~np.isnan(M).any(axis=1) & np.isfinite(case)
    usable &= (np.abs(case) < args.clip) & (np.abs(M) < args.clip).all(axis=1)
    if usable.sum() < 1000:
        sys.exit(f"Only {usable.sum()} bins usable across the panel and case")

    Mu = np.clip(M[usable], -args.clip, args.clip)
    cu = np.clip(case[usable], -args.clip, args.clip)

    # Centre each panel sample, then each bin: the per-bin step is what a
    # pooled reference already does, and the decomposition then sees only
    # what that leaves behind.
    Mu = Mu - np.median(Mu, axis=0, keepdims=True)
    bin_med = np.median(Mu, axis=1)
    R = Mu - bin_med[:, None]

    # Left singular vectors span the bin-space directions along which
    # panel samples vary together -- the correlated bias.
    U, S, _ = np.linalg.svd(R, full_matrices=False)
    k = min(args.components, U.shape[1])
    Uk = U[:, :k]

    r = cu - np.median(cu) - bin_med
    denoised_resid = r - Uk @ (Uk.T @ r)

    out = np.array(case, dtype=float)
    out[usable] = denoised_resid
    out[~usable] = np.nan

    kept = 0
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for row, val in zip(rows, out):
            if not np.isfinite(val):
                # A bin the panel could not inform is dropped rather than
                # passed through uncorrected next to corrected neighbours.
                continue
            row = dict(row)
            row["log2"] = f"{val:.6f}"
            w.writerow(row)
            kept += 1
    var_removed = 1 - np.var(denoised_resid) / np.var(r) if np.var(r) else 0.0
    print(f"{os.path.basename(args.case)}: panel={len(panel_paths)} k={k} "
          f"bins={kept}/{len(rows)} variance_removed={var_removed:.1%}")


if __name__ == "__main__":
    main()
