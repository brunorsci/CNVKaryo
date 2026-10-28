rule gatk_contig_ploidy:
    """
    Per-chromosome ploidy for every sample. This is the whole-chromosome
    arm of the workflow: a trisomy appears here, not in the focal caller,
    and the ploidy estimate is also what GermlineCNVCaller conditions on
    so that a focal call is made against the right baseline.

    The priors table states, per contig, the prior probability of each
    ploidy state. It is intentionally permissive on the autosomes so that
    a genuine trisomy is reachable rather than penalised to zero.
    """
    input:
        filtered = config["output"] + "gcnv/filtered.interval_list",
        counts = expand(config["output"] + "gcnv/counts/{sample}.counts.hdf5", sample=SAMPLES),
        priors = config["gcnv"]["ploidy_priors"]
    output:
        calls = directory(config["output"] + "gcnv/ploidy-calls"),
        model = directory(config["output"] + "gcnv/ploidy-model")
    params:
        counts_args = lambda wc, input: " ".join(f"-I {c}" for c in input.counts),
        prefix = config["output"] + "gcnv/ploidy"
    conda:
        "../envs/gatk.yaml"
    log:
        config["output"] + "logs/gcnv/contig_ploidy.log"
    shell:
        """
        gatk DetermineGermlineContigPloidy \
            -L {input.filtered} --interval-merging-rule OVERLAPPING_ONLY \
            {params.counts_args} \
            --contig-ploidy-priors {input.priors} \
            --output $(dirname {params.prefix}) \
            --output-prefix ploidy 2>&1 | tee {log}
        """


rule gatk_germline_cnv_caller:
    """
    Cohort-mode focal CNV calling. Bias removal and copy-number inference
    happen inside one model, which is the design claim of the method: the
    authors argue that separate normalisation and calling steps can
    discard real CNV signal.
    """
    input:
        filtered = config["output"] + "gcnv/filtered.interval_list",
        annotated = config["output"] + "gcnv/annotated_intervals.tsv",
        counts = expand(config["output"] + "gcnv/counts/{sample}.counts.hdf5", sample=SAMPLES),
        ploidy = config["output"] + "gcnv/ploidy-calls"
    output:
        calls = directory(config["output"] + "gcnv/cohort-calls"),
        model = directory(config["output"] + "gcnv/cohort-model")
    params:
        counts_args = lambda wc, input: " ".join(f"-I {c}" for c in input.counts),
        outdir = config["output"] + "gcnv"
    threads: 16
    conda:
        "../envs/gatk.yaml"
    log:
        config["output"] + "logs/gcnv/germline_cnv_caller.log"
    shell:
        """
        gatk GermlineCNVCaller \
            --run-mode COHORT \
            -L {input.filtered} --interval-merging-rule OVERLAPPING_ONLY \
            {params.counts_args} \
            --contig-ploidy-calls {input.ploidy} \
            --annotated-intervals {input.annotated} \
            --output {params.outdir} \
            --output-prefix cohort 2>&1 | tee {log}
        """


rule gatk_postprocess_calls:
    """
    Per-sample VCFs plus a per-interval copy-number track. The segments
    VCF is what the karyogram and the evidence tables are built from.
    """
    input:
        calls = config["output"] + "gcnv/cohort-calls",
        model = config["output"] + "gcnv/cohort-model",
        ploidy = config["output"] + "gcnv/ploidy-calls"
    output:
        intervals_vcf = config["output"] + "gcnv/vcf/{sample}_intervals.vcf.gz",
        segments_vcf = config["output"] + "gcnv/vcf/{sample}_segments.vcf.gz",
        ratios = config["output"] + "gcnv/vcf/{sample}_denoised_copy_ratios.tsv"
    params:
        index = lambda wc: SAMPLES.index(wc.sample)
    conda:
        "../envs/gatk.yaml"
    log:
        config["output"] + "logs/gcnv/{sample}_postprocess.log"
    shell:
        """
        gatk PostprocessGermlineCNVCalls \
            --calls-shard-path {input.calls} \
            --model-shard-path {input.model} \
            --contig-ploidy-calls {input.ploidy} \
            --sample-index {params.index} \
            --output-genotyped-intervals {output.intervals_vcf} \
            --output-genotyped-segments {output.segments_vcf} \
            --output-denoised-copy-ratios {output.ratios} \
            2>&1 | tee {log}
        """


rule gatk_gcnv_case_mode:
    """
    Call a sample against a model trained earlier, instead of training a
    new one.

    This is what makes the method usable on one sample at a time. COHORT
    mode needs at least 30 exomes (GATK's documented minimum; ~200 is the
    suggested batch size, not a requirement), but the model it produces
    is reusable: CASE mode has no lower limit on sample count.

    The model is an asset per capture kit, not per run. GATK requires the
    case sample to be "strictly compatible (in terms of library
    preparation and sequencing protocol)" with the training cohort, and
    the same capture kit for exomes. A sample from a kit with no trained
    model cannot be called this way -- either 30 samples of that kit are
    assembled and a model trained, or a different caller is used. The
    pipeline fails loudly rather than calling against a mismatched model.

    Note that supplying a model in COHORT mode does not reuse it: GATK
    uses it only to initialise and then re-parametrises from the input
    counts. Reuse means CASE mode.
    """
    input:
        counts = config["output"] + "gcnv/counts/{sample}.counts.hdf5",
        model = lambda wc: config["gcnv"]["pretrained_model"],
        ploidy = config["output"] + "gcnv/ploidy-calls"
    output:
        calls = directory(config["output"] + "gcnv/case-calls/{sample}")
    params:
        outdir = config["output"] + "gcnv/case-calls"
    threads: 8
    log:
        config["output"] + "logs/gcnv/{sample}_case_mode.log"
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        gatk GermlineCNVCaller \
            --run-mode CASE \
            -I {input.counts} \
            --contig-ploidy-calls {input.ploidy} \
            --model {input.model} \
            --output {params.outdir} \
            --output-prefix {wildcards.sample} 2>&1 | tee {log}
        """


rule gatk_export_model:
    """
    Copies the trained cohort model out of the run directory into a
    named, versioned location so later runs can call against it in CASE
    mode. Without this the model is buried in a results tree that gets
    deleted with the rest of a run -- which has already happened once in
    this project, to a set of alignments that were not wrong, only the
    results built on them.
    """
    input:
        model = config["output"] + "gcnv/cohort-model",
        ploidy_model = config["output"] + "gcnv/ploidy-model"
    output:
        stamp = config["output"] + "gcnv/exported_model/MODEL_INFO.txt"
    params:
        dest = config["output"] + "gcnv/exported_model",
        kit = config.get("gcnv", {}).get("kit_name", "unspecified_kit"),
        intervals = config["output"] + "gcnv/filtered.interval_list"
    log:
        config["output"] + "logs/gcnv/export_model.log"
    shell:
        """
        mkdir -p {params.dest}
        cp -r {input.model} {params.dest}/cohort-model
        cp -r {input.ploidy_model} {params.dest}/ploidy-model
        cp {params.intervals} {params.dest}/filtered.interval_list
        {{
          echo "capture_kit: {params.kit}"
          echo "trained: $(date -Iseconds)"
          echo "gatk: $(gatk --version 2>&1 | grep -oE 'v?4\\.[0-9.]+' | head -1)"
          echo "n_training_samples: $(ls {input.model}/.. 2>/dev/null | wc -l)"
          echo "NOTE: a case sample must share library prep, sequencing protocol"
          echo "      and capture kit with the cohort this model was trained on."
        }} > {output.stamp}
        cat {output.stamp} | tee {log}
        """
