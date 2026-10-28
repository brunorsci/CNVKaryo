rule cnvkit_access:
    """
    Genome-level accessible-region BED for WGS coverage binning (excludes
    assembly gaps). Sample-independent -- built once and reused by every
    WGS sample in the run.
    """
    input:
        fasta = config["genome_assets"][config["genome_build"]]["ref"]
    output:
        access_bed = config["output"] + "output_cnvkaryo/cnvkit/access.bed"
    log:
        config["output"] + "logs/cnvkit_reference/cnvkit_access.log"
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        cnvkit.py access {input.fasta} -o {output.access_bed} 2>&1 | tee {log}
        """
