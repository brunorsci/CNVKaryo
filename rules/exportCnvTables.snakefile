rule export_cnv_tables:
    """
    Produces the three audit tables for the virtual karyotype
    (segments/bands/breakpoints) from {sample}.call.cns -- see
    scripts/export_cnv_tables.py. {sample}_cnv_bands.tsv becomes the single
    drawing source for every plot (ideogram, circular, genome-wide and
    per-chromosome scatter): one band row renders as exactly one colored
    bar in the image.
    """
    input:
        call_cns=config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.call.cns",
        sex=config["output"] + "output_cnvkaryo/cnvkit/{sample}_sex.txt",
        cytoband=config["genome_assets"][config["genome_build"]]["cytoband"]
    output:
        segments=config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_segments.tsv",
        bands=config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_bands.tsv",
        breakpoints=config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_breakpoints.tsv"
    params:
        declared_sex=lambda wc: config.get("sample_sex", {}).get(wc.sample, "")
    log:
        config["output"] + "logs/export_cnv_tables/{sample}.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        DECLARED="{params.declared_sex}"
        python3 scripts/export_cnv_tables.py \
            ${{DECLARED:+--declared-sex $DECLARED}} \
            --call-cns {input.call_cns} \
            --sex-file {input.sex} \
            --cytoband {input.cytoband} \
            --sample {wildcards.sample} \
            --out-segments {output.segments} \
            --out-bands {output.bands} \
            --out-breakpoints {output.breakpoints} \
            2>&1 | tee {log}
        """
