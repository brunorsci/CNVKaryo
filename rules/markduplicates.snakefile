rule markDuplicates:
    """
    Marks AND removes PCR/optical duplicates (REMOVE_DUPLICATES=true) --
    required for depth-based CNV calling: duplicates inflate coverage
    non-uniformly (correlated with GC content and repeat structure), which
    would otherwise bias the log2 ratio CNVkit computes per bin.

    Feeds CNVkit directly: cnvkaryo has no variant-calling step, so there
    is no BQSR here (unlike the original VKaryo). BaseRecalibrator only
    rewrites per-base QUAL scores -- it does not change read position,
    CIGAR, MAPQ or duplicate status, so it cannot affect a depth-based
    tool's output. Running it would only add a dbSNP dependency and
    15-30 minutes per sample for zero effect on the CNV calls.
    """
    input:
        config["output"] + "output_cnvkaryo/aligned/{sample}_aligned.bam",
    output:
        bam=config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam",
        bai=config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam.bai",
        txt=config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates_metrics.txt",
    conda:
        "../envs/picard.yaml"
    resources:
        mem_mb=25000,
    benchmark: config["output"] + "output_cnvkaryo/benchmarks/mark_duplicates/marked_{sample}_benchmark.txt",
    threads: 4
    shell:
        """
        picard MarkDuplicates \
        I={input} \
        O={output.bam}\
        M={output.txt} \
        REMOVE_DUPLICATES=true \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ADD_PG_TAG_TO_READS=false \
        PROGRAM_RECORD_ID=null \
        ASSUME_SORTED=true

        picard BuildBamIndex INPUT={output.bam} OUTPUT={output.bai}
        """
