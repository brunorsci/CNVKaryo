rule fastp:
    input:
        r1=config["input"]["fastq_dir"] + "{sample}_1.fastq.gz",
        r2=config["input"]["fastq_dir"] + "{sample}_2.fastq.gz"
    output:
        r1=config["output"] + "output_cnvkaryo/fastp/{sample}_fastp_1.fastq.gz",
        r2=config["output"] + "output_cnvkaryo/fastp/{sample}_fastp_2.fastq.gz",
        html=config["output"] + "output_cnvkaryo/fastp/{sample}_fastp.html",
        json=config["output"] + "output_cnvkaryo/fastp/{sample}_fastp.json",
    log:
        config["output"] + "output_cnvkaryo/fastp/log/{sample}.log"
    threads: 6
    conda:
        "../envs/fastp.yaml"
    benchmark: config["output"] + "output_cnvkaryo/benchmarks/fastp/fastp_{sample}_benchmark.txt",
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --html {output.html} --json {output.json}
        """
