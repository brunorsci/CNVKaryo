rule cnvkit_reference:
    """
    Builds the flat reference for hybrid-capture (WES) coverage, with
    GC/RepeatMasker correction. Do not add --no-gc/--no-rmask: without that
    correction a flat reference produces a systematic genome-wide log2
    shift per sample, which silently breaks every downstream CNV call.
    """
    input:
        target_bed=config["output"] + "output_cnvkaryo/cnvkit/capture_kit.target.bed",
        antitarget_bed=config["output"] + "output_cnvkaryo/cnvkit/capture_kit.antitarget.bed",
        fasta=config["genome_assets"][config["genome_build"]]["ref"]
    output:
        reference=config["output"] + "output_cnvkaryo/cnvkit/flat_reference.cnn"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_reference.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py reference \
            -t {input.target_bed} \
            -a {input.antitarget_bed} \
            -f {input.fasta} \
            -o {output.reference} \
            2>&1 | tee {log}
        """
