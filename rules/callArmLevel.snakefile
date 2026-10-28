rule call_arm_level:
    """
    Whole-chromosome and arm-level calls read straight from the per-bin
    log2 ratios, which is what a conventional karyogram reports. This is
    deliberately independent of the CBS segments the band tables are built
    from: segmentation optimises breakpoint placement, so on noisy data it
    emits short high-variance segments and calling copy number per segment
    reports an aneuploidy wherever noise crosses a threshold. See
    scripts/call_arm_level.py for the magnitude + direction-consistency
    criteria.

    Sex chromosomes are excluded unless --call-sex-chroms is passed; a
    documented sex from config["sample_sex"] overrides cnvkit.py sex, for
    the same reason as in cnvkit_call.
    """
    input:
        cnr = config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cnr",
        sex = config["output"] + "output_cnvkaryo/cnvkit/{sample}_sex.txt",
        cytoband = config["genome_assets"][config["genome_build"]]["cytoband"],
        # Declared so that editing the calling logic re-triggers the call.
        # Cheap rule: it reads one .cnr and writes two small text files.
        script = "scripts/call_arm_level.py",
        common = "scripts/cnvkaryo_plot_common.py"
    output:
        arms = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_arm_level.tsv",
        karyotype = config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_karyotype.txt"
    params:
        sample = "{sample}",
        declared_sex = lambda wc: config.get("sample_sex", {}).get(wc.sample, "")
    log:
        config["output"] + "logs/call_arm_level/{sample}.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        DECLARED=""
        if [ -n "{params.declared_sex}" ]; then
            DECLARED="--declared-sex {params.declared_sex}"
        fi
        python {input.script} \
            --cnr {input.cnr} \
            --cytoband {input.cytoband} \
            --sex-file {input.sex} \
            $DECLARED \
            --sample {params.sample} \
            --out-arms {output.arms} \
            --out-karyotype {output.karyotype} \
            2>&1 | tee {log}
        """
