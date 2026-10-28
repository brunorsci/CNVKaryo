rule cnv_scatter_plot:
    """
    Per-chromosome CNV scatter (log2 per bin) for the chromosome carrying a
    sample's gene of interest, with the band colored by state (gain/loss/
    neutral, from absolute cn in {sample}_cnv_bands.tsv) and the gene label
    anchored on the band that actually overlaps it -- see
    scripts/plot_cnv_scatter_colored.py. Only runs for genes explicitly
    listed in config["cnv"]["genes_of_interest"] -- never for an arbitrary
    gene by default (see scripts/plot_gene.sh for on-demand single-gene
    requests, which register the gene in that list first).
    """
    input:
        cnr = config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cnr",
        bands = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_bands.tsv",
        gene_coords = config["genome_assets"][config["genome_build"]]["gene_coords"]
    output:
        png = config["output"] + "output_cnvkaryo/karyogram/{sample}/genes/{sample}_{gene}_scatter.png"
    params:
        chrom = lambda wc: GENE_TO_CHROM[(wc.sample, wc.gene)],
        title = lambda wc: f"{GENE_TO_CELLLINE.get((wc.sample, wc.gene), wc.sample)} ({wc.sample}) -- {GENE_TO_CHROM[(wc.sample, wc.gene)]}",
    log:
        config["output"] + "logs/cnv_scatter/{sample}_{gene}.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        python3 scripts/plot_cnv_scatter_colored.py \
            --cnr {input.cnr} \
            --bands {input.bands} \
            --chrom {params.chrom} \
            --gene {wildcards.gene} \
            --gene-coords {input.gene_coords} \
            --output {output.png} \
            --title "{params.title}" \
            2>&1 | tee {log}
        """


def _genome_scatter_gene_args(wc):
    genes = SAMPLE_TO_GENES.get(wc.sample, [])
    return " ".join(f"--gene {g}" for g in genes)


CHROM_LIST_ALL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


rule cnv_genome_scatter_plot:
    """
    Virtual karyotype: genome-wide 24-chromosome grid plus 24 individual
    high-resolution per-chromosome images, in one pass -- see
    scripts/plot_cnv_genome_scatter.py. Color/classification source is
    {sample}_cnv_bands.tsv (absolute cn, sex-aware): one band row always
    renders as exactly one colored bar. The karyotype call from
    {sample}_arm_level.tsv is shaded behind them, so each panel shows the
    verdict and the segment evidence it was read from at once. Runs for every sample in this run;
    any configured gene(s) of interest for that sample are highlighted too.
    """
    input:
        cnr = config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cnr",
        bands = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_bands.tsv",
        arms = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_arm_level.tsv",
        gene_coords = config["genome_assets"][config["genome_build"]]["gene_coords"],
        script = "scripts/plot_cnv_genome_scatter.py",
        common = "scripts/cnvkaryo_plot_common.py"
    output:
        png = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_genome_wide.png",
        chromosomes = expand(
            config["output"] + "output_cnvkaryo/karyogram/{{sample}}/chromosomes/{chrom}.png",
            chrom=CHROM_LIST_ALL
        )
    params:
        title = lambda wc: f"{SAMPLE_TO_CELLLINE.get(wc.sample, wc.sample)} ({wc.sample}) -- virtual karyotype",
        gene_args = _genome_scatter_gene_args,
        chrom_dir = config["output"] + "output_cnvkaryo/karyogram/{sample}/chromosomes"
    log:
        config["output"] + "logs/cnv_scatter/{sample}_genome_wide.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        python3 {input.script} \
            --bands {input.bands} \
            --arms {input.arms} \
            --cnr {input.cnr} \
            --gene-coords {input.gene_coords} \
            {params.gene_args} \
            --output {output.png} \
            --output-chromosomes-dir {params.chrom_dir} \
            --title "{params.title}" \
            2>&1 | tee {log}
        """
