# Reference files

cnvkaryo does not ship reference data. Three things are needed.

## 1. Genome (GRCh38) and BWA index

```bash
mkdir -p reference && cd reference
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

The BWA index takes a few hours and writes alongside the FASTA.

Only needed when starting from FASTQ. If you already have aligned BAMs, point
the config at them and skip this.

## 2. Capture target BED

For exome data, the target BED of the kit that prepared the libraries. Kit
vendors distribute these; for whole-genome input, a BED of the chromosomes is
enough and the workflow can generate one from the FASTA index.

The expected-coverage statistics are kit-specific: a BED from a different kit
produces wrong expectations, not merely noisier ones.

## 3. Cytobands

Used to place calls on the ideogram, in the standard UCSC format:

```bash
wget -O reference/cytoBand.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
gunzip reference/cytoBand.txt.gz
```

## Checking

Every path is set in `config/config.yaml`. To confirm they all resolve without
running anything:

```bash
snakemake --configfile config/example.yaml --rulegraph > /dev/null && echo OK
```
