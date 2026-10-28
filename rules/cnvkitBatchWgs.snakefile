rule cnvkit_batch_wgs:
    """
    WGS coverage, flat reference and segmentation in a single pass. Unlike
    hybrid-capture mode, WGS needs no target/antitarget BED pair: bins are
    derived directly from the accessible-genome BED, and `-n` with no
    filenames tells CNVkit to build a flat reference on the fly rather than
    reuse a pre-built pooled reference. Produces the same {sample}.cnr/.cns
    naming as the hybrid-mode rule, so every downstream rule (cnvkit_sex,
    cnvkit_call, export_cnv_tables, the plots) is agnostic to which
    sequencing type produced them.
    """
    input:
        bam=config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam",
        bai=config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam.bai",
        fasta=config["genome_assets"][config["genome_build"]]["ref"],
        access_bed=config["output"] + "output_cnvkaryo/cnvkit/access.bed",
        refflat=config["genome_assets"][config["genome_build"]]["refflat"]
    output:
        cnr=config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cnr",
        cns=config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.cns"
    params:
        output_dir=config["output"] + "output_cnvkaryo/cnvkit/"
    log:
        config["output"] + "logs/cnvkit_batch/{sample}.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 10
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

        cnvkit.py batch {input.bam} \
            -n \
            -m wgs \
            -f {input.fasta} \
            -g {input.access_bed} \
            --annotate {input.refflat} \
            --output-dir {params.output_dir} \
            --processes {threads} \
            2>&1 | tee {log}
        """
