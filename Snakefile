"""
cnvkaryo - Virtual Karyotyping from WES/WGS CNV calls
========================================================

Standalone pipeline extracted from VKaryo (virtual karyotype module):
fastp -> bwa-mem -> MarkDuplicates -> CNVkit -> cnvkit.py call (absolute
copy number, sex-aware) -> audit tables (segments/bands/breakpoints) ->
plots (ideogram, circular plot, aneuploidy summary, per-gene scatter,
per-chromosome and genome-wide virtual karyotype). No RNA-seq dependency.

No BQSR: cnvkaryo does no variant calling, and BaseRecalibrator only
rewrites per-base QUAL scores, which a depth-based CNV caller never reads
-- it would add a dbSNP dependency and real runtime for zero effect on
the CNV calls.

Every plot's color/classification comes from the absolute copy number in
{sample}.call.cns via {sample}_cnv_bands.tsv, never an independent log2
cutoff per script.
"""

import os
import sys
import csv

configfile: "config.yaml"

SAMPLES = config.get("samples")
if not SAMPLES:
    (SAMPLES,) = glob_wildcards(config["input"]["fastq_dir"] + "{sample}_1.fastq.gz")

print(f"cnvkaryo: {len(SAMPLES)} sample(s): {SAMPLES}")

SEQUENCING_TYPE = config.get("sequencing_type", "wes")
if SEQUENCING_TYPE not in ("wes", "wgs"):
    raise ValueError(f"Unknown sequencing_type {SEQUENCING_TYPE!r}: must be 'wes' or 'wgs'")

# Which CNVkit reference a WES sample is called against. A flat reference's
# bin geometry is baked in from whatever capture kit built it -- reusing
# one cohort's kit-specific reference (e.g. the Agilent SureSelect V6
# LL-100 panel) for a sample captured with a DIFFERENT kit produces
# biologically impossible calls (cn in the dozens, log2 swinging +/-20),
# because coverage gets counted over bins that don't match where that
# sample's probes actually sit. Three reference modes, picked per sample:
#   - "agilent_v6": the original LL-100-panel flat reference (default --
#     exact previous behavior, for samples actually captured with that kit)
#   - "generic_exome": a kit-agnostic flat reference built from the union
#     of RefSeq exons (rules/cnvkitGenericReference.snakefile) -- for
#     standalone external samples whose real capture kit is unknown
#   - paired: for a tumor with a configured matched normal
#     (cnv_reference.paired_normals), the reference is built from that
#     normal's own real coverage instead of any flat model (standard
#     CNVkit tumor/normal practice -- see cnvkitPairedReference.snakefile)
_cnv_reference_cfg = config.get("cnv_reference", {})
_cnv_reference_default_mode = _cnv_reference_cfg.get("default_mode", "agilent_v6")
TUMOR_TO_NORMAL = dict(_cnv_reference_cfg.get("paired_normals", {}))
NORMAL_SAMPLES = set(TUMOR_TO_NORMAL.values())


def cnvkit_reference_for(wildcards):
    sample = wildcards.sample
    if sample in TUMOR_TO_NORMAL:
        normal = TUMOR_TO_NORMAL[sample]
        return config["output"] + f"output_cnvkaryo/cnvkit/{normal}_paired_reference.cnn"
    if sample in NORMAL_SAMPLES:
        return config["output"] + f"output_cnvkaryo/cnvkit/{sample}_paired_reference.cnn"
    if _cnv_reference_default_mode == "generic_exome":
        return config["output"] + "output_cnvkaryo/cnvkit/flat_reference_generic.cnn"
    if _cnv_reference_default_mode == "agilent_v6":
        return config["output"] + "output_cnvkaryo/cnvkit/flat_reference.cnn"
    if _cnv_reference_default_mode == "pooled":
        if sample in set(_cnv_reference_cfg.get("panel_samples", [])):
            raise ValueError(
                f"{sample} is both a test sample and a member of the pooled panel. "
                f"A sample corrected by a panel containing itself has part of its own "
                f"copy number absorbed into the expected value.")
        return config["output"] + "output_cnvkaryo/cnvkit/pooled_reference.cnn"
    raise ValueError(f"Unknown cnv_reference.default_mode {_cnv_reference_default_mode!r}")


# Which caller performs copy-number inference. GATK-gCNV is wrapped, not
# reimplemented: it carries the best published performance in the size
# regime of the recurrent microdeletions (0.99 recall / 0.96 precision
# for deletions spanning >=10 exons). It needs a cohort of kit-matched
# samples, which an ordinary single-sample run does not have, so it is
# selected per run rather than always on.
CNV_CALLER = config.get("cnv_caller", "cnvkit")
if CNV_CALLER not in ("cnvkit", "gatk_gcnv", "both"):
    raise ValueError(f"cnv_caller must be cnvkit, gatk_gcnv or both, not {CNV_CALLER!r}")
if CNV_CALLER in ("gatk_gcnv", "both"):
    if SEQUENCING_TYPE != "wes":
        raise ValueError("gatk_gcnv is a hybrid-capture workflow; use cnvkit for wgs")
    if "gcnv" not in config or "ploidy_priors" not in config["gcnv"]:
        raise ValueError("gatk_gcnv needs config['gcnv']['ploidy_priors']")


include: "rules/fastp.snakefile"
include: "rules/bwa_mem.snakefile"
include: "rules/markduplicates.snakefile"

# Coverage/reference-building differs between hybrid-capture (WES) and WGS;
# everything downstream of {sample}.cnr/.cns is identical either way.
_run_cnvkit = CNV_CALLER in ("cnvkit", "both")
if SEQUENCING_TYPE == "wes" and _run_cnvkit:
    if _cnv_reference_default_mode == "agilent_v6":
        include: "rules/cnvkitTargets.snakefile"
        include: "rules/cnvkitRef.snakefile"
    if _cnv_reference_default_mode == "generic_exome" or TUMOR_TO_NORMAL:
        include: "rules/cnvkitGenericReference.snakefile"
    if _cnv_reference_default_mode == "pooled":
        include: "rules/cnvkitTargets.snakefile"
        include: "rules/cnvkitPooledReference.snakefile"
    if TUMOR_TO_NORMAL:
        include: "rules/cnvkitPairedReference.snakefile"
    include: "rules/cnvkitBatch.snakefile"
elif _run_cnvkit:
    include: "rules/cnvkitAccess.snakefile"
    include: "rules/cnvkitBatchWgs.snakefile"

include: "rules/cnvkitCall.snakefile"
include: "rules/exportCnvTables.snakefile"
include: "rules/cnv_scatter.snakefile"
include: "rules/plot_aneuploidy.snakefile"
include: "rules/callArmLevel.snakefile"

# Needs a panel to tell a biallelic deletion from a region the capture
# never reaches, so it only runs where one is configured.
if _cnv_reference_default_mode == "pooled":
    include: "rules/callHomozygousDeletions.snakefile"

# The karyogram is drawn from whichever ratio track a run produced. Both
# are plotted when both callers run, so the figure can be compared
# across callers on identical axes.
KARYOGRAM_SOURCES = []
if CNV_CALLER in ("cnvkit", "both"):
    KARYOGRAM_SOURCES.append("cnvkit")
if CNV_CALLER in ("gatk_gcnv", "both"):
    KARYOGRAM_SOURCES.append("gcnv")


def _karyogram_ratio_input(wildcards):
    if wildcards.source == "cnvkit":
        return config["output"] + f"output_cnvkaryo/cnvkit/{wildcards.sample}_markduplicates.cnr"
    return config["output"] + f"gcnv/vcf/{wildcards.sample}_denoised_copy_ratios.tsv"


include: "rules/karyogram.snakefile"


if CNV_CALLER in ("gatk_gcnv", "both"):
    include: "rules/gatkGcnvIntervals.snakefile"
    include: "rules/gatkGcnvCall.snakefile"

    # cohort: train a model from >=30 kit-matched samples and export it.
    # case:   call against a model trained earlier, any number of samples.
    GCNV_MODE = config.get("gcnv", {}).get("mode", "cohort")
    if GCNV_MODE not in ("cohort", "case"):
        raise ValueError(f"gcnv.mode must be cohort or case, not {GCNV_MODE!r}")
    if GCNV_MODE == "cohort" and len(SAMPLES) < 30:
        print(f"WARNING: gCNV cohort mode with {len(SAMPLES)} samples. GATK documents "
              f"a minimum of 30 for exomes; the bias model will be poorly determined "
              f"below that and the calls should not be trusted.", file=sys.stderr)
    if GCNV_MODE == "case" and not config.get("gcnv", {}).get("pretrained_model"):
        raise ValueError("gcnv.mode=case needs gcnv.pretrained_model pointing at a "
                         "model trained on the SAME capture kit")


if config.get("validation", {}).get("truth_dir"):
    include: "rules/evaluateValidation.snakefile"

# Per-gene scatter plots only run for genes explicitly listed in
# config["cnv"]["genes_of_interest"] -- never "any gene in the genome" by
# default. Restricted to this run's SAMPLES.
GENE_TO_CHROM = {}
GENE_TO_CELLLINE = {}
SAMPLE_TO_GENES = {}
SAMPLE_TO_CELLLINE = {}
CNV_SCATTER_OUTPUTS = []
_genes_of_interest_path = config.get("cnv", {}).get("genes_of_interest")
if _genes_of_interest_path and os.path.exists(_genes_of_interest_path):
    with open(_genes_of_interest_path) as _f:
        for _row in csv.DictReader(_f, delimiter="\t"):
            if _row["dna_sample"] not in SAMPLES:
                continue
            _key = (_row["dna_sample"], _row["gene"])
            GENE_TO_CHROM[_key] = _row["chrom"]
            GENE_TO_CELLLINE[_key] = _row["cell_line"]
            SAMPLE_TO_GENES.setdefault(_row["dna_sample"], []).append(_row["gene"])
            SAMPLE_TO_CELLLINE[_row["dna_sample"]] = _row["cell_line"]
            CNV_SCATTER_OUTPUTS.append(
                config["output"] + f"output_cnvkaryo/karyogram/{_row['dna_sample']}/genes/{_row['dna_sample']}_{_row['gene']}_scatter.png"
            )
else:
    print(f"WARNING: {_genes_of_interest_path} not found -- skipping per-gene-of-interest scatter plots.")

# Virtual karyotype genome-wide grid (runs for every sample, with or without
# a documented gene of interest) -- the 24 per-chromosome images are
# declared outputs of the same rule (cnv_genome_scatter_plot).
CNV_GENOME_SCATTER_OUTPUTS = [
    config["output"] + f"output_cnvkaryo/karyogram/{s}/{s}_genome_wide.png" for s in SAMPLES
]

KARYOGRAM_OUTPUTS = (
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_aneuploidy_ideogram.pdf", sample=SAMPLES) +
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_circular.pdf", sample=SAMPLES) +
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_aneuploidy_summary.txt", sample=SAMPLES) +
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_segments.tsv", sample=SAMPLES) +
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_bands.tsv", sample=SAMPLES) +
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_cnv_breakpoints.tsv", sample=SAMPLES) +
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_arm_level.tsv", sample=SAMPLES) +
    expand(config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_karyotype.txt", sample=SAMPLES) +
    CNV_SCATTER_OUTPUTS +
    CNV_GENOME_SCATTER_OUTPUTS
)

if _cnv_reference_default_mode == "pooled":
    KARYOGRAM_OUTPUTS = KARYOGRAM_OUTPUTS + expand(
        config["output"] + "output_cnvkaryo/karyogram/{sample}/{sample}_homozygous_deletions.tsv",
        sample=SAMPLES)

KARYOGRAM_OUTPUTS = KARYOGRAM_OUTPUTS + expand(
    config["output"] + "karyogram/{sample}/{sample}_genome_trace_{source}.png",
    sample=SAMPLES, source=KARYOGRAM_SOURCES)

if CNV_CALLER in ("gatk_gcnv", "both") and GCNV_MODE == "cohort":
    KARYOGRAM_OUTPUTS = KARYOGRAM_OUTPUTS + [
        config["output"] + "gcnv/exported_model/MODEL_INFO.txt"]
if CNV_CALLER in ("gatk_gcnv", "both"):
    KARYOGRAM_OUTPUTS = KARYOGRAM_OUTPUTS + expand(
        config["output"] + "gcnv/vcf/{sample}_segments.vcf.gz", sample=SAMPLES) + expand(
        config["output"] + "gcnv/vcf/{sample}_denoised_copy_ratios.tsv", sample=SAMPLES) + [
        config["output"] + "gcnv/ploidy-calls"]

# Scoring against documented karyotypes only applies to a run that has
# them; an ordinary diagnostic run has nothing to score against.
if config.get("validation", {}).get("truth_dir"):
    KARYOGRAM_OUTPUTS = KARYOGRAM_OUTPUTS + [
        config["output"] + "output_cnvkaryo/validation/validation_report.txt"]

rule all:
    input:
        KARYOGRAM_OUTPUTS
