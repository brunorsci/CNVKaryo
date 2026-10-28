rule plot_aneuploidy:
    """
    Schematic ideogram + circular plot + aneuploidy text summary, driven by
    {sample}_arm_level.tsv. The CBS bands are drawn underneath as a
    detail layer but no longer decide the verdict: a chromosome-scale
    judgement read off segment lengths reported 156 non-neutral bands on
    a constitutional trisomy 21 whose only true abnormality is +21.
    """
    input:
        bands = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_bands.tsv",
        arms = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_arm_level.tsv",
        script = "scripts/plot_aneuploidy_ideogram.py"
    output:
        ideogram = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_aneuploidy_ideogram.pdf",
        circular = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_circular.pdf",
        summary = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_aneuploidy_summary.txt"
    params:
        sample = "{sample}"
    log:
        config["output"] + "logs/plot_aneuploidy/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        python {input.script} \
            --bands {input.bands} \
            --arms {input.arms} \
            --sample {params.sample} \
            --output {output.ideogram} \
            --output-circular {output.circular} \
            --output-summary {output.summary} \
            2>&1 | tee {log}
        """
