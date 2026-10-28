rule gatk_preprocess_intervals:
    """
    Capture intervals prepared for GATK-gCNV. For exome data the bins are
    the capture targets themselves, so binning is switched off
    (--bin-length 0) and the targets are padded, which recovers reads
    falling just outside the probe footprint.
    """
    input:
        bed = config["genome_assets"][config["genome_build"]]["exome_bed_noalt"],
        ref = config["genome_assets"][config["genome_build"]]["ref"],
        dict_ = config["genome_assets"][config["genome_build"]]["ref"].replace(".fasta", ".dict")
    output:
        intervals = config["output"] + "gcnv/preprocessed.interval_list"
    params:
        padding = config.get("gcnv", {}).get("padding", 250)
    conda:
        "../envs/gatk.yaml"
    log:
        config["output"] + "logs/gcnv/preprocess_intervals.log"
    shell:
        """
        gatk PreprocessIntervals \
            -L {input.bed} -R {input.ref} \
            --bin-length 0 --padding {params.padding} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {output.intervals} 2>&1 | tee {log}
        """


rule gatk_annotate_intervals:
    """
    GC content per interval. The GATK team recommend this step even though
    it is optional: without it the cohort model has no GC covariate, and
    GC is the dominant systematic bias in hybrid capture.
    """
    input:
        intervals = config["output"] + "gcnv/preprocessed.interval_list",
        ref = config["genome_assets"][config["genome_build"]]["ref"]
    output:
        annotated = config["output"] + "gcnv/annotated_intervals.tsv"
    conda:
        "../envs/gatk.yaml"
    log:
        config["output"] + "logs/gcnv/annotate_intervals.log"
    shell:
        """
        gatk AnnotateIntervals \
            -L {input.intervals} -R {input.ref} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {output.annotated} 2>&1 | tee {log}
        """


rule gatk_collect_read_counts:
    """
    Per-interval read counts, the only sample-level input gCNV consumes.
    HDF5 so the cohort step can read many samples cheaply.
    """
    input:
        bam = config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam",
        bai = config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam.bai",
        intervals = config["output"] + "gcnv/preprocessed.interval_list",
        ref = config["genome_assets"][config["genome_build"]]["ref"]
    output:
        counts = config["output"] + "gcnv/counts/{sample}.counts.hdf5"
    conda:
        "../envs/gatk.yaml"
    log:
        config["output"] + "logs/gcnv/{sample}_counts.log"
    shell:
        """
        gatk CollectReadCounts \
            -I {input.bam} -L {input.intervals} -R {input.ref} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --format HDF5 -O {output.counts} 2>&1 | tee {log}
        """


rule gatk_filter_intervals:
    """
    Drops intervals that the cohort cannot inform: extreme GC, low
    mappability, and low or highly variable coverage across samples. The
    gCNV paper filters intervals whose median count falls below ten reads
    per sample within a batch, for the same reason this pipeline drops
    zero-depth bins -- an interval with no reads carries no evidence and
    only contributes noise.
    """
    input:
        intervals = config["output"] + "gcnv/preprocessed.interval_list",
        annotated = config["output"] + "gcnv/annotated_intervals.tsv",
        counts = expand(config["output"] + "gcnv/counts/{sample}.counts.hdf5", sample=SAMPLES)
    output:
        filtered = config["output"] + "gcnv/filtered.interval_list"
    params:
        counts_args = lambda wc, input: " ".join(f"-I {c}" for c in input.counts)
    conda:
        "../envs/gatk.yaml"
    log:
        config["output"] + "logs/gcnv/filter_intervals.log"
    shell:
        """
        gatk FilterIntervals \
            -L {input.intervals} --annotated-intervals {input.annotated} \
            {params.counts_args} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {output.filtered} 2>&1 | tee {log}
        """
