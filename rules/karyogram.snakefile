rule genome_trace:
    """
    The virtual karyogram: one continuous genomic axis, every informative
    interval as a faint point, and the median over fixed genomic windows
    as a line.

    Nothing on the figure comes from a call. An aneuploidy or a large
    deletion is visible because the point cloud moves, which is the only
    thing that makes the figure evidence rather than an illustration of a
    conclusion already reached. Calls live in the tables beside it.

    Also writes one panel per chromosome. The genome-wide view compresses
    a chromosome into a few dozen pixels, which shows a whole-chromosome
    shift and hides an arm-level or sub-arm one; the per-chromosome
    panels use a finer window and mark the centromere, because a shift
    that stops at the centromere is a different finding from one that
    crosses it.

    Drawn from whichever ratio track the run produced -- CNVkit .cnr or
    the GATK denoised copy ratios -- so the two callers can be compared
    on identical axes.
    """
    input:
        ratios = _karyogram_ratio_input,
        script = "scripts/plot_genome_trace.py",
        common = "scripts/cnvkaryo_plot_common.py",
        cytoband = config["genome_assets"][config["genome_build"]]["cytoband"]
    output:
        png = config["output"] + "karyogram/{sample}/{sample}_genome_trace_{source}.png",
        chroms = directory(config["output"] + "karyogram/{sample}/chromosomes_{source}")
    params:
        window_mb = config.get("karyogram", {}).get("window_mb", 1.0),
        ylim = config.get("karyogram", {}).get("ylim", 2.0),
        chrom_window_mb = config.get("karyogram", {}).get("chrom_window_mb", 0.25),
        title = lambda wc: f"{wc.sample} ({wc.source})"
    log:
        config["output"] + "logs/karyogram/{sample}_{source}_trace.log"
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        python {input.script} \
            --ratios {input.ratios} \
            --sample {wildcards.sample} \
            --title "{params.title}" \
            --cytoband {input.cytoband} \
            --window-mb {params.window_mb} \
            --ylim {params.ylim} \
            --per-chromosome-dir {output.chroms} \
            --per-chromosome-window-mb {params.chrom_window_mb} \
            --out {output.png} \
            2>&1 | tee {log}
        """
