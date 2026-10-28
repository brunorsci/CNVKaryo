rule filter_alt_contig_probes:
    """
    Drops capture-kit probes referencing alternate/random/unplaced contigs
    that the pipeline's reference FASTA does not carry (see
    scripts/filter_alt_contig_probes.py) -- extracting sequence for those
    contigs during GC/RepeatMasker correction would otherwise raise a
    KeyError.
    """
    input:
        bed = config["genome_assets"][config["genome_build"]]["exome_bed"],
        fai = config["genome_assets"][config["genome_build"]]["ref"] + ".fai"
    output:
        noalt_bed = config["output"] + "output_cnvkaryo/cnvkit/capture_kit.noalt.bed"
    log:
        config["output"] + "logs/cnvkit_reference/filter_alt_contig_probes.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/filter_alt_contig_probes.py \
            --bed {input.bed} \
            --fai {input.fai} \
            --output {output.noalt_bed} \
            2>&1 | tee {log}
        """


rule cnvkit_target:
    input:
        bed = config["output"] + "output_cnvkaryo/cnvkit/capture_kit.noalt.bed",
        refflat = config["genome_assets"][config["genome_build"]]["refflat"]
    output:
        target_bed = config["output"] + "output_cnvkaryo/cnvkit/capture_kit.target.bed"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_target.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py target {input.bed} \
            --annotate {input.refflat} \
            --split \
            -o {output.target_bed} \
            2>&1 | tee {log}
        """


rule cnvkit_antitarget:
    input:
        target_bed = config["output"] + "output_cnvkaryo/cnvkit/capture_kit.target.bed"
    output:
        antitarget_bed = config["output"] + "output_cnvkaryo/cnvkit/capture_kit.antitarget.bed"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_antitarget.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py antitarget {input.target_bed} \
            -o {output.antitarget_bed} \
            2>&1 | tee {log}
        """
