#!/usr/bin/env python3
"""
Drops capture-kit BED probes referencing alternate/random/unplaced contigs
("_alt", "_random", "chrUn_...") that the pipeline's reference FASTA does
not carry.

Without this filter, GC/RepeatMasker correction in the cnvkit_reference
rule raises a KeyError on the missing contigs, which would otherwise force
running with `--no-gc --no-rmask` (an uncorrected flat reference) -- and
that produces a systematic genome-wide log2 shift per sample, silently
breaking every downstream CNV call.

Usage:
    python3 filter_alt_contig_probes.py \
        --bed reference/capture_kit.bed \
        --fai reference/genome.fasta.fai \
        --output reference/capture_kit.noalt.bed
"""

import argparse


def load_valid_contigs(fai_path):
    contigs = set()
    with open(fai_path) as f:
        for line in f:
            contigs.add(line.split('\t')[0])
    return contigs


def filter_bed(bed_path, valid_contigs, output_path):
    kept, dropped = 0, 0
    with open(bed_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#'):
                fout.write(line)
                continue
            chrom = line.split('\t', 1)[0]
            if chrom in valid_contigs:
                fout.write(line)
                kept += 1
            else:
                dropped += 1
    return kept, dropped


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bed', required=True)
    ap.add_argument('--fai', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    valid = load_valid_contigs(args.fai)
    kept, dropped = filter_bed(args.bed, valid, args.output)
    print(f"Valid contigs in reference FASTA: {len(valid)}")
    print(f"Probes kept: {kept}")
    print(f"Probes dropped (contig missing from FASTA): {dropped}")
