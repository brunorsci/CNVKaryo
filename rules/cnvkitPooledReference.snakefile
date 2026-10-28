rule cnvkit_pooled_reference:
    """
    Reference built from the observed coverage of a named panel of
    samples, rather than from a flat model.

    CNVkit's documentation is explicit that a pooled reference supplies,
    per bin, an expected value AND a spread measured across the panel,
    and that the inverse squared spread is what weights each bin during
    segmentation. A flat reference carries no spread, so neither the
    weighting nor the systematic-bias removal happens; the documentation
    states that with a flat reference, exome noise "is expected". That
    was measured here: window noise fell from 0.141 to 0.101 on the same
    samples and bins when only the reference changed.

    Panel membership comes from config["cnv_reference"]["panel_samples"]
    and is a deliberate choice, not every available sample. The panel's
    per-bin median is contaminated by whatever copy-number changes its
    members actually carry, so members are selected for being as close
    to diploid and as structurally simple as the cohort allows. A sample
    under test must never appear in the panel that corrects it.
    """
    input:
        targets = lambda wc: expand(
            config["output"] + "output_cnvkaryo/cnvkit/{s}_markduplicates.targetcoverage.cnn",
            s=config["cnv_reference"]["panel_samples"]),
        antitargets = lambda wc: expand(
            config["output"] + "output_cnvkaryo/cnvkit/{s}_markduplicates.antitargetcoverage.cnn",
            s=config["cnv_reference"]["panel_samples"]),
        fasta = config["genome_assets"][config["genome_build"]]["ref"]
    output:
        reference = config["output"] + "output_cnvkaryo/cnvkit/pooled_reference.cnn"
    log:
        config["output"] + "logs/cnvkit_reference/pooled_reference.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py reference {input.targets} {input.antitargets} \
            -f {input.fasta} -o {output.reference} 2>&1 | tee {log}
        """


rule cnvkit_panel_coverage:
    """
    Coverage for a panel member.

    Panel samples are not analysed, so no other rule computes their
    coverage: `cnvkit batch` only runs for samples under test. This rule
    exists so the panel can be assembled from samples that are never
    themselves called, which is what keeps a sample out of the reference
    that corrects it.
    """
    input:
        bam = config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam",
        bai = config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam.bai",
        target = config["output"] + "output_cnvkaryo/cnvkit/capture_kit.target.bed",
        antitarget = config["output"] + "output_cnvkaryo/cnvkit/capture_kit.antitarget.bed"
    output:
        target_cnn = config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.targetcoverage.cnn",
        antitarget_cnn = config["output"] + "output_cnvkaryo/cnvkit/{sample}_markduplicates.antitargetcoverage.cnn"
    log:
        config["output"] + "logs/cnvkit_reference/{sample}_panel_coverage.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 4
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py coverage {input.bam} {input.target} -p {threads} -o {output.target_cnn} 2>&1 | tee {log}
        cnvkit.py coverage {input.bam} {input.antitarget} -p {threads} -o {output.antitarget_cnn} 2>&1 | tee -a {log}
        """
