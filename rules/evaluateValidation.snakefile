rule evaluate_validation:
    """
    Scores the called karyotypes against the documented ones and writes the
    confusion matrix with PPV/NPV and Wilson intervals.

    Only included when config["validation"]["truth_dir"] is set, because it
    needs documented karyotypes -- it is a property of the validation
    panel, not of an ordinary run. It lives in the workflow rather than
    being computed by hand so the published numbers are reproducible from
    the repository.
    """
    input:
        arms = expand(config["output"] +
                      "output_cnvkaryo/karyogram/{sample}/{sample}_arm_level.tsv",
                      sample=SAMPLES),
        script = "scripts/evaluate_karyotype.py",
        cytoband = config["genome_assets"][config["genome_build"]]["cytoband"]
    output:
        report = config["output"] + "output_cnvkaryo/validation/validation_report.txt"
    params:
        karyogram_dir = config["output"] + "output_cnvkaryo/karyogram",
        truth_dir = config["validation"]["truth_dir"],
        strata = config["validation"].get("strata", ""),
        min_coverage = config["validation"].get("min_coverage", 0.5)
    log:
        config["output"] + "logs/evaluate_validation/evaluate.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        STRATA=""
        if [ -n "{params.strata}" ]; then STRATA="--strata {params.strata}"; fi
        python {input.script} \
            --karyogram-dir {params.karyogram_dir} \
            --truth-dir {params.truth_dir} \
            --cytoband {input.cytoband} \
            --min-coverage {params.min_coverage} \
            $STRATA \
            --out {output.report} \
            2>&1 | tee {log}
        """
