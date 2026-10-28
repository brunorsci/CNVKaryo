# Running the example

The quickest way to confirm the installation works is to run cnvkaryo on a
public sample with a known answer. No data ships with this repository.

## The sample

**SRR1190472** — whole-exome sequencing of an individual with **trisomy 21**,
from BioProject `PRJNA240916` (Chaiyasap et al., whole-genome and exome
sequencing of monozygotic twins discordant for trisomy 21).

It is a good first run because the expected result is unambiguous: chromosome
21 should come out at three copies and every other chromosome at two. If the
ideogram shows that, the installation is correct.

## Downloading

With the SRA Toolkit:

```bash
mkdir -p example/data && cd example/data
prefetch SRR1190472
fasterq-dump --split-files SRR1190472
gzip SRR1190472_1.fastq SRR1190472_2.fastq
```

Or directly from ENA, which serves the FASTQ without the toolkit:

```bash
mkdir -p example/data && cd example/data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/002/SRR1190472/SRR1190472_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/002/SRR1190472/SRR1190472_2.fastq.gz
```

About 6 GB. Rename to the convention the workflow expects:

```bash
mv SRR1190472_1.fastq.gz sample1_1.fastq.gz
mv SRR1190472_2.fastq.gz sample1_2.fastq.gz
```

## Reference files

cnvkaryo needs the GRCh38 genome with a BWA index, and the capture target BED
of the kit used. See `docs/references.md`.

## Running

```bash
cp config/config.yaml config/example.yaml
# edit: input directory, output directory, reference paths

snakemake --configfile config/example.yaml --cores 8 --use-conda
```

Runtime depends on the machine; on 8 cores expect a few hours, most of it in
alignment.

## What to look for

| Output | Expected |
|---|---|
| `results/ideogram/sample1_aneuploidy.png` | chromosome 21 marked as gain, everything else neutral |
| `results/karyotype/sample1_genome.png` | the chr21 trace sitting one step above the rest |
| `results/tables/sample1_chromosome_calls.tsv` | `chr21` called `gain`, ratio near 1.5 |

A ratio near **1.5** is three copies against a two-copy baseline, which is what
trisomy looks like in read depth.

## Other public samples

`docs/validation.md` lists further public accessions used during development,
including a matched tumour/normal pair and a complex-karyotype AML case, if you
want to exercise the pipeline on harder material.
