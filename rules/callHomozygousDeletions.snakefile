rule call_homozygous_deletions:
    """
    Biallelic losses, detected as runs of empty intervals where the
    reference panel has coverage.

    This is the one event class the rest of the pipeline is structurally
    blind to. A homozygous deletion drives coverage to zero, and the log2
    ratio of zero coverage is meaningless -- CNVkit emits about -24 --
    so every log2-based step either discards the interval or is poisoned
    by it. The arm-level caller drops zero-depth intervals explicitly.
    The signal is therefore read from raw depth, against the panel.
    """
    input:
        coverage = config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.targetcoverage.cnn",
        panel = lambda wc: expand(
            config["output"] + "output_cnvkaryo/cnvkit/{s}_markduplicates.targetcoverage.cnn",
            s=config["cnv_reference"]["panel_samples"]),
        ratios = config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cnr",
        script = "scripts/call_homozygous_deletions.py"
    output:
        calls = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_homozygous_deletions.tsv"
    log:
        config["output"] + "logs/homozygous_deletions/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        python {input.script} \
            --coverage {input.coverage} \
            --ratios {input.ratios} \
            --panel-coverage {input.panel} \
            --sample {wildcards.sample} \
            --out {output.calls} \
            2>&1 | tee {log}
        """
