rule generic_exon_bed:
    """
    Kit-agnostic target BED (union of RefSeq exons) for samples whose real
    capture kit is unknown -- see scripts/refflat_to_exon_bed.py.
    """
    input:
        refflat=config["genome_assets"][config["genome_build"]]["refflat"]
    output:
        bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.bed"
    log:
        config["output"] + "logs/cnvkit_reference/generic_exon_bed.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python3 scripts/refflat_to_exon_bed.py \
            --refflat {input.refflat} \
            --output {output.bed} \
            2>&1 | tee {log}
        """


rule cnvkit_generic_target:
    input:
        bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.bed",
        refflat=config["genome_assets"][config["genome_build"]]["refflat"]
    output:
        target_bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.target.bed"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_generic_target.log"
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


rule cnvkit_generic_antitarget:
    input:
        target_bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.target.bed"
    output:
        antitarget_bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.antitarget.bed"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_generic_antitarget.log"
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


rule cnvkit_generic_reference:
    """
    Flat reference over the kit-agnostic exon BED -- used for samples with
    no matched normal and an unknown/non-Agilent capture kit (see
    cnvkit_reference_for() in the Snakefile).
    """
    input:
        target_bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.target.bed",
        antitarget_bed=config["output"] + "output_cnvkaryo/cnvkit/generic_exons.antitarget.bed",
        fasta=config["genome_assets"][config["genome_build"]]["ref"]
    output:
        reference=config["output"] + "output_cnvkaryo/cnvkit/flat_reference_generic.cnn"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_generic_reference.log"
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
