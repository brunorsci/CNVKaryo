rule cnvkit_batch:
    input:
        bam=config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam",
        bai=config["output"] + "output_cnvkaryo/marked/{sample}_markduplicates.bam.bai",
        reference=cnvkit_reference_for
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
        # Force conda's libstdc++ over the system's (system lacks CXXABI_1.3.15)
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

        cnvkit.py batch {input.bam} \
            -r {input.reference} \
            --output-dir {params.output_dir} \
            --processes {threads} \
            --method hybrid \
            --drop-low-coverage \
            2>&1 | tee {log}
        """
