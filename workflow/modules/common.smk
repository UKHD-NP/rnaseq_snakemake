TRUE_VALUES = {"true", "yes", "1", "t", "y"}


def as_bool(value, default=False):
    """
    Normalize config values to boolean.
    
    Accepts various representations of 'True':
    - True (boolean)
    - 'true' (lowercase string)
    - 'True' (capitalized string)
    - 'TRUE' (uppercase string)
    - 1 (integer)
    Returns False for other values.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in TRUE_VALUES
    if isinstance(value, int):
        return value == 1
    return default


def is_enabled(module_name, default=False):
    """
    Check if a module is enabled in the config.

    - Section missing or empty -> returns *default*.
    - Section present but no 'enabled' key -> True.
    - Section present with 'enabled' key -> as_bool(value).

    Examples:
        is_enabled("trimming")                   # off when section missing
        is_enabled("bam_filter", default=True)    # on unless explicitly disabled
    """
    module_cfg = config.get(module_name, {})
    if not module_cfg:
        return default
    if isinstance(module_cfg, dict):
        if "enabled" not in module_cfg:
            return True
        return as_bool(module_cfg["enabled"], default=default)
    return as_bool(module_cfg, default=default)


def get_sample_rows(sample_id):
    """Return all samplesheet rows matching a sample ID."""
    sample_rows = samplesheet[samplesheet['sample_id'] == sample_id]
    if sample_rows.empty:
        fatal(f"Sample ID '{sample_id}' not found in samplesheet")
    return sample_rows


def get_outdir(sample_id):
    """Get the output directory for a given sample ID."""
    sample_rows = get_sample_rows(sample_id)
    outdirs = sample_rows['outdir'].dropna().unique().tolist()
    if len(outdirs) != 1:
        fatal(
            f"Sample ID '{sample_id}' maps to multiple outdir values in samplesheet: {outdirs}. "
            "Please keep outdir consistent for duplicated sample_id rows."
        )
    return outdirs[0]


def get_raw_lane_fastqs(wildcards):
    """
    Return all raw FASTQ lane files for a sample ID.
    Supports multiple rows per sample_id in the samplesheet.
    """
    sample_rows = get_sample_rows(wildcards.sample_id)
    fq1s = sample_rows['fq1'].tolist()
    fq2s = sample_rows['fq2'].tolist()
    if len(fq1s) != len(fq2s):
        fatal(f"Sample ID '{wildcards.sample_id}' has unequal fq1/fq2 counts: {len(fq1s)} vs {len(fq2s)}")

    invalid_fq1 = [str(path) for path in fq1s if not str(path).strip().lower().endswith(".gz")]
    invalid_fq2 = [str(path) for path in fq2s if not str(path).strip().lower().endswith(".gz")]
    if invalid_fq1 or invalid_fq2:
        invalid_parts = []
        if invalid_fq1:
            invalid_parts.append(f"fq1={invalid_fq1}")
        if invalid_fq2:
            invalid_parts.append(f"fq2={invalid_fq2}")
        fatal(
            f"Sample ID '{wildcards.sample_id}' requires gzipped FASTQ inputs (*.gz); "
            f"found non-gz paths: {'; '.join(invalid_parts)}"
        )

    return fq1s, fq2s


def get_raw_lane_fq1(wildcards):
    return get_raw_lane_fastqs(wildcards)[0]


def get_raw_lane_fq2(wildcards):
    return get_raw_lane_fastqs(wildcards)[1]


# Merge raw FASTQ lanes for duplicated sample IDs.
# For samples with one lane, this is a simple copy-through via cat.
rule merge_raw_fastqs:
    input:
        fq1 = get_raw_lane_fq1,
        fq2 = get_raw_lane_fq2
    output:
        fq1 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_2.fastq.gz")
    log:
        os.path.join("{outdir}", "logs", "raw_merge", "{sample_id}.merge_raw_fastqs.log")
    message:
        "{wildcards.sample_id}: Merging raw FASTQ lanes"
    shell:
        """
        mkdir -p $(dirname {output.fq1})
        mkdir -p $(dirname {log})
        rm -f {output.fq1} {output.fq2}

        N_R1=$(echo {input.fq1} | wc -w)
        N_R2=$(echo {input.fq2} | wc -w)

        if [ "$N_R1" -eq 1 ] && [ "$N_R2" -eq 1 ]; then
            # Single-lane sample: use symlink to avoid duplicate storage
            ln -sf "$(readlink -f {input.fq1})" {output.fq1}
            ln -sf "$(readlink -f {input.fq2})" {output.fq2}
            echo "[INFO] Mode: symlink (single lane)" > {log}
        else
            # Multi-lane sample: concatenate lanes in listed order
            cat {input.fq1} > {output.fq1}
            cat {input.fq2} > {output.fq2}
            echo "[INFO] Mode: merge (multi-lane)" > {log}
        fi

        if [ ! -s "{output.fq1}" ] || [ ! -s "{output.fq2}" ]; then
            echo "[ERROR] Merged FASTQ output is empty." >> {log}
            exit 1
        fi

        echo "[INFO] Merged R1 inputs: {input.fq1}" >> {log}
        echo "[INFO] Merged R2 inputs: {input.fq2}" >> {log}
        """


def get_paired_fq(wildcards):
    """Get merged raw FASTQ paths for a sample."""
    outdir = get_outdir(wildcards.sample_id)
    return [
        os.path.join(outdir, "raw_merged", f"{wildcards.sample_id}_merged_{read}.fastq.gz")
        for read in ("1", "2")
    ]


def get_paired_trimmed_fq(wildcards):
    """Return trimmed FASTQ paths if trimming is enabled, otherwise raw FASTQs."""
    if is_enabled("trimming"):
        # Check if outdir is available in wildcards, otherwise resolve from samplesheet
        outdir = getattr(wildcards, 'outdir', get_outdir(wildcards.sample_id))

        # Both fastp and trim_galore output to the same "trim" directory
        trim_dir = os.path.join(outdir, "trim")
        return [
            os.path.join(trim_dir, f"{wildcards.sample_id}_trimmed_{read}.fastq.gz")
            for read in ("1", "2")
        ]
    else:
        return get_paired_fq(wildcards)


def get_paired_sortmerna_fq(wildcards):
    """Return SortMeRNA non-rRNA FASTQ paths for a sample."""
    outdir = getattr(wildcards, 'outdir', get_outdir(wildcards.sample_id))
    sortmerna_dir = os.path.join(outdir, "sortmerna")
    return [
        os.path.join(sortmerna_dir, f"{wildcards.sample_id}.nonrna_{read}.fastq.gz")
        for read in ("1", "2")
    ]


def get_paired_align_fq(wildcards):
    """Return FASTQ paths feeding alignment/fusion: SortMeRNA output if enabled, else trimmed/raw FASTQs."""
    if is_enabled("sortmerna"):
        return get_paired_sortmerna_fq(wildcards)
    return get_paired_trimmed_fq(wildcards)


def get_bam_basename(sample_id):
    """Return BAM basename based on markduplicates setting."""
    if is_enabled("markduplicates"):
        return f"{sample_id}.markdup.sorted.bam"
    return f"{sample_id}.bam"


def get_bam(wildcards):
    """Return BAM path chosen for downstream rules."""
    return os.path.join(wildcards.outdir, "bam", get_bam_basename(wildcards.sample_id))


def get_bam_bai(wildcards):
    """Return the BAM index path matching get_bam()."""
    return f"{get_bam(wildcards)}.bai"


def get_ref_bed():
    """Return BED reference path used by RSeQC rules."""
    return config.get("ref", {}).get("bed", "")


def is_rseqc_submodule_enabled(submodule_name):
    """Check whether a specific RSeQC submodule is enabled."""
    return (
        is_enabled("rseqc")
        and config.get("rseqc", {}).get(submodule_name, {}).get("enabled", False)
    )


def get_rseqc_targets(outdir, sample_id, purpose="workflow"):
    """
    Build RSeQC target paths for either workflow completion or MultiQC inputs.
    """
    if purpose not in {"workflow", "multiqc"}:
        fatal(f"Unsupported RSeQC target purpose: {purpose}")

    junction_annotation = (
        os.path.join(outdir, "logs", "rseqc", f"{sample_id}.junction_annotation.log")
        if purpose == "multiqc"
        else os.path.join(outdir, "rseqc", "junction_annotation", f"{sample_id}.junction.bed")
    )

    submodule_targets = {
        "bam_stat": os.path.join(outdir, "rseqc", "bam_stat", f"{sample_id}.bam_stat.txt"),
        "infer_experiment": os.path.join(outdir, "rseqc", "infer_experiment", f"{sample_id}.infer_experiment.txt"),
        "read_distribution": os.path.join(outdir, "rseqc", "read_distribution", f"{sample_id}.read_distribution.txt"),
        "inner_distance": os.path.join(outdir, "rseqc", "inner_distance", f"{sample_id}.inner_distance.txt"),
        "junction_annotation": junction_annotation,
        "tin": os.path.join(outdir, "rseqc", "tin", f"{sample_id}.tin.summary.txt"),
    }

    if purpose == "workflow":
        submodule_targets.update({
            "gene_body_coverage": os.path.join(outdir, "rseqc", "gene_body_coverage", f"{sample_id}.geneBodyCoverage.curves.pdf"),
            "read_GC": os.path.join(outdir, "rseqc", "read_GC", f"{sample_id}.GC_plot.pdf"),
            "read_duplication": os.path.join(outdir, "rseqc", "read_duplication", f"{sample_id}.DupRate_plot.pdf"),
            "junction_saturation": os.path.join(outdir, "rseqc", "junction_saturation", f"{sample_id}.junctionSaturation_plot.pdf"),
        })

    return [
        path for submodule_name, path in submodule_targets.items()
        if is_rseqc_submodule_enabled(submodule_name)
    ]


def get_target_files(sample_ids):
    """Determine target files for the workflow based on enabled modules."""
    targets = []
    
    # Process each sample
    for sample_id in sample_ids:
        outdir = get_outdir(sample_id)

        def _path(*parts):
            return os.path.join(outdir, *parts)
        
        # Always include MultiQC report
        targets.append(_path("multiqc", f"{sample_id}.multiqc.html"))

        # Module-specific outputs
        conditional_targets = (
            ("salmon_counts", _path("salmon", f"{sample_id}.quant.sf")),
            ("featurecounts", _path("featurecounts", f"{sample_id}.fc")),
            ("fusion", _path("arriba", f"{sample_id}.fusions.tsv")),
            ("stringtie", _path("stringtie", f"{sample_id}.ballgown")),
            ("dupradar", _path("dupradar", f"{sample_id}_duprateExpDens.pdf")),
        )
        for module_name, path in conditional_targets:
            if is_enabled(module_name):
                targets.append(path)

        if is_enabled("bigwig"):
            targets.extend([
                _path("bigwig", f"{sample_id}.forward.CPM.bw"),
                _path("bigwig", f"{sample_id}.reverse.CPM.bw"),
            ])

        if is_enabled("rseqc"):
            targets.extend(get_rseqc_targets(outdir, sample_id, purpose="workflow"))
            
        if is_enabled("genome_load_keep_memory"):
            targets.append(_path("bam", f"{sample_id}.star_memory_removal.log"))

        # Always include deletion log
        targets.append(_path("logs", f"{sample_id}.deletion.log"))
    
    return targets
