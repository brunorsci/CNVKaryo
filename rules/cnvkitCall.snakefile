rule cnvkit_sex:
    """
    Infers chromosomal sex from the sample's own coverage (cnvkit.py sex).
    Only used when config["sample_sex"] has no entry for the sample: the
    heuristic compares chrX/chrY coverage, so a tumour that lost one X
    reads as male (observed on HCC1395, a documented female donor, where
    chrX LOH gave X log2=-0.95 against Y log2=-20, i.e. no Y reads at all).
    A documented sex always beats this inference -- see cnvkit_call.
    """
    input:
        cnr=config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cnr"
    output:
        sex=config["output"] + "output_cnvkaryo/cnvkit/{sample}_sex.txt"
    log:
        config["output"] + "logs/cnvkit_call/{sample}_sex.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        cnvkit.py sex {input.cnr} -o {output.sex} 2>&1 | tee {log}
        """


rule cnvkit_call:
    """
    Explicit absolute copy number call (cnvkit.py call -x <sex>) rather than
    cnvkit.py batch's implicit default -- {sample}.call.cns becomes the
    single source of truth for every downstream classification/color (see
    scripts/export_cnv_tables.py and cnvkaryo_plot_common.py), not the raw
    .cns, which carries an uncalibrated per-chromosome log2 offset.
    """
    input:
        cns=config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cns",
        sex=config["output"] + "output_cnvkaryo/cnvkit/{sample}_sex.txt"
    output:
        call_cns=config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.call.cns"
    params:
        declared_sex=lambda wc: config.get("sample_sex", {}).get(wc.sample, "")
    log:
        config["output"] + "logs/cnvkit_call/{sample}_call.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        if [ -n "{params.declared_sex}" ]; then
            SAMPLE_SEX="{params.declared_sex}"
            echo "Declared sex (config sample_sex): $SAMPLE_SEX" | tee {log}
        else
            SAMPLE_SEX=$(tail -n1 {input.sex} | cut -f2)
            echo "Inferred sex (cnvkit.py sex): $SAMPLE_SEX" | tee {log}
        fi
        # --filter ci merges adjacent segments that are statistically
        # indistinguishable (overlapping confidence interval) before
        # calling copy number -- without it, CBS over-segmentation produces
        # many spurious focal gain/loss calls.
        # --center median re-centers log2 on the median before applying the
        # absolute cn cutoffs -- without it, a per-sample global log2 offset
        # (e.g. a slight GC/depth mismatch against the flat reference) can
        # push nearly the whole genome to the wrong side of the cutoff
        # (log2=-0.4 everywhere => cn=1 everywhere => a false genome-wide
        # monosomy call). cnvkit.py batch applies both internally by
        # default; they are made explicit here because this rule now calls
        # `call` on its own rather than relying on batch's implicit default.
        # No -y here: the references this pipeline builds already carry a
        # diploid-chrX baseline, so telling `call` to assume a male
        # reference makes it read a neutral chrX (log2 ~ 0) as a single
        # copy. Measured directly on the self-referenced HCC1395BL normal,
        # where chrX cannot be anything but neutral by construction:
        #   -x Female     -> chrX cn=2  (correct)
        #   -x Female -y  -> chrX cn=1  (spurious monosomy)
        cnvkit.py call {input.cns} --filter ci --center median --drop-low-coverage -x "$SAMPLE_SEX" -o {output.call_cns} 2>&1 | tee -a {log}
        """
