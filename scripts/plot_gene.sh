#!/usr/bin/env bash
# Requests the CNV scatter plot for ONE specific gene in one sample -- the
# official mechanism for investigating a gene outside the automatic
# gene-of-interest list. Not a loose bypass: it registers the (sample, gene)
# pair in the configured list (config["cnv"]["genes_of_interest"]) before
# asking Snakemake for the plot, so why that gene was looked at stays on
# record.
#
# Usage: scripts/plot_gene.sh <sample> <GENE> [cell_line] [expected_alteration]
set -euo pipefail

CNVKARYO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$CNVKARYO_DIR"

SAMPLE="${1:?Usage: scripts/plot_gene.sh <sample> <GENE> [cell_line] [expected_alteration]}"
GENE="${2:?Usage: scripts/plot_gene.sh <sample> <GENE> [cell_line] [expected_alteration]}"
CELL_LINE="${3:-$SAMPLE}"
EXPECTED="${4:-manual_investigation}"

GENE_COORDS="reference/gene_coords_hg38.tsv"
GENES_OF_INTEREST="reference/genes_of_interest.tsv"

CHROM=$(awk -F'\t' -v g="$GENE" '$1==g {print $2; exit}' "$GENE_COORDS" || true)
if [ -z "$CHROM" ]; then
    echo "ERROR: gene '$GENE' not found in $GENE_COORDS -- check the symbol" >&2
    exit 1
fi

if ! awk -F'\t' -v s="$SAMPLE" -v g="$GENE" 'NR>1 && $2==s && $3==g {found=1} END{exit !found}' "$GENES_OF_INTEREST"; then
    printf '%s\t%s\t%s\t%s\t%s\n' "$CELL_LINE" "$SAMPLE" "$GENE" "$CHROM" "$EXPECTED" >> "$GENES_OF_INTEREST"
    echo "Registered in $GENES_OF_INTEREST: sample=$SAMPLE gene=$GENE chrom=$CHROM"
else
    echo "Pair sample=$SAMPLE gene=$GENE was already registered in $GENES_OF_INTEREST"
fi

OUTPUT_DIR=$(python3 -c "import yaml; c=yaml.safe_load(open('config.yaml')); print(c['output'])")
TARGET="${OUTPUT_DIR}output_cnvkaryo/karyogram/${SAMPLE}/genes/${SAMPLE}_${GENE}_scatter.png"

echo "Requesting: $TARGET"
snakemake --cores 4 -p "$TARGET"
