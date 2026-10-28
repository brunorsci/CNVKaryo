rule cnvkit_paired_reference:
    """
    Reference built from a matched normal's own coverage rather than a
    flat/generic model -- the standard CNVkit tumor/normal workflow. Used
    for pairs configured in cnv_reference.paired_normals (Snakefile), so a
    tumor's actual capture kit never needs to be known: the normal's own
    coverage under whatever kit was really used defines the baseline.

    Still needs a target/antitarget BED to fix bin boundaries, since none
    is sample-specific here; reuses the kit-agnostic exon BED (see
    cnvkitGenericReference.snakefile) for that purpose only -- the real
    bias correction comes from the paired normal's coverage, not from the
    BED matching any particular kit's bait design.

    The configured normal MUST be the same chromosomal sex as the sample
    it references, even for a same-batch pseudo-pair that isn't a true
    clinical tumor/normal from one individual: a cross-sex reference's
    chrX coverage is a large outlier that drags the genome-wide median
    used for centering (cnvkit_call's --center median) with it, producing
    a uniform false gain/loss across every autosome -- a result that
    looks clean (few, large segments) and is actually wrong, which is
    worse than the flat reference's obvious noise.
    """
    input:
        normal_bam=config["output"] + "output_cnvkaryo/marked/{normal_sample}_markduplicates.bam",
        normal_bai=config["output"] + "output_cnvkaryo/marked/{normal_sample}_markduplicates.bam.bai",
        target_bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.target.bed",
        antitarget_bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.antitarget.bed",
        fasta=config["genome_assets"][config["genome_build"]]["ref"],
        refflat=config["genome_assets"][config["genome_build"]]["refflat"]
    output:
        reference=config["output"] + "output_cnvkaryo/cnvkit/{normal_sample}_paired_reference.cnn"
    params:
        output_dir=config["output"] + "output_cnvkaryo/cnvkit/"
    log:
        config["output"] + "logs/cnvkit_reference/{normal_sample}_paired_reference.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 10
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

        cnvkit.py batch \
            -n {input.normal_bam} \
            -f {input.fasta} \
            -t {input.target_bed} \
            -a {input.antitarget_bed} \
            --annotate {input.refflat} \
            --output-reference {output.reference} \
            --output-dir {params.output_dir} \
            --processes {threads} \
            2>&1 | tee {log}
        """
