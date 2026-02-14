TRUE_VALUES = {"true", "yes", "1", "t", "y"}


def error_msg(msg):
    """Format parser/config errors consistently across modules."""
    return f"[ERROR] {msg}"


def as_bool(value, default=False):
    """
    Normalize config values to boolean.
    Accepts bool, common truthy strings, and integer 1.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in TRUE_VALUES
    if isinstance(value, int):
        return value == 1
    return default


def is_enabled(module_name):
    """
    Check if a module is enabled in the config.
    Accepts various representations of 'True':
    - True (boolean)
    - 'true' (lowercase string)
    - 'True' (capitalized string)
    - 'TRUE' (uppercase string)
    - 1 (integer)
    Returns False for other values.
    """
    # Check if module exists in config with an 'enabled' key
    if module_name not in config or "enabled" not in config[module_name]:
        return False
    
    return as_bool(config[module_name]["enabled"], default=False)

def get_sample_rows(sample_id):
    """Return all samplesheet rows matching a sample ID."""
    sample_rows = samplesheet[samplesheet['sample_id'] == sample_id]
    if sample_rows.empty:
        raise ValueError(error_msg(f"Sample ID '{sample_id}' not found in samplesheet"))
    return sample_rows

# Get the output directory for a given sample id
def get_outdir(sample_id):
    """Get the output directory for a given sample ID with error handling."""
    sample_rows = get_sample_rows(sample_id)
    outdirs = sample_rows['outdir'].dropna().unique().tolist()
    if len(outdirs) != 1:
        raise ValueError(
            error_msg(
                f"Sample ID '{sample_id}' maps to multiple outdir values in samplesheet: {outdirs}. "
                "Please keep outdir consistent for duplicated sample_id rows."
            )
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
        raise ValueError(
            error_msg(
                f"Sample ID '{wildcards.sample_id}' has unequal fq1/fq2 counts: {len(fq1s)} vs {len(fq2s)}"
            )
        )
    return fq1s, fq2s

def get_raw_lane_fq1(wildcards):
    return get_raw_lane_fastqs(wildcards)[0]

def get_raw_lane_fq2(wildcards):
    return get_raw_lane_fastqs(wildcards)[1]

# Get list of paths to raw fastq files from samplesheet
def get_paired_fq(wildcards):
    """Get merged raw FASTQ paths for a sample."""
    outdir = get_outdir(wildcards.sample_id)
    return [
        os.path.join(outdir, "raw_merged", f"{wildcards.sample_id}_merged_{read}.fastq.gz")
        for read in ("1", "2")
    ]

# Get list of paths to trimmed fastq files (if trimming is disabled, return raw fastq files)
def get_paired_trimmed_fq(wildcards):
    """
    Returns paths to trimmed FASTQ files if trimming is enabled,
    otherwise returns paths to the original raw FASTQ files.
    Supports both fastp and trim_galore.
    """
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

def get_ref_bed():
    """Return BED reference path used by RSeQC rules."""
    return config.get("ref", {}).get("bed", "")


def get_bam_basename(sample_id):
    """Return BAM basename based on markduplicates setting."""
    if is_enabled("markduplicates"):
        return f"{sample_id}.markdup.sorted.bam"
    return f"{sample_id}.bam"


def get_bam(wildcards):
    """Return BAM path chosen for downstream rules."""
    return os.path.join(wildcards.outdir, "bam", get_bam_basename(wildcards.sample_id))

def get_bam_bai(wildcards):
    """
    Returns the BAM index path matching get_bam().
    """
    return f"{get_bam(wildcards)}.bai"

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
        raise ValueError(error_msg(f"Unsupported RSeQC target purpose: {purpose}"))

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
        
# Collect MultiQC inputs based on which modules are enabled in the config
def get_input_multiqc(wildcards):
    """
    Collect input files for MultiQC analysis based on configuration settings.
    
    This function dynamically builds a list of input files for MultiQC, 
    including core files and conditionally adding module-specific outputs.
    
    Args:
        wildcards: Snakemake wildcards object containing sample information
    
    Returns:
        list: Paths to files to be processed by MultiQC
    """
    sample_id = wildcards.sample_id
    outdir = get_outdir(sample_id)

    def _path(*parts):
        return os.path.join(outdir, *parts)

    # Core mandatory outputs
    targets = [_path("bam", f"{sample_id}.Log.final.out")]

    # Conditionally add trimming quality control reports
    if is_enabled("trimming"):
        trimming_tool = config.get('trimming', {}).get('tool', 'fastp')

        # Add raw FastQC reports (before trimming) - for both tools
        targets.extend([
            _path("fastqc_raw", f"{sample_id}_raw_1_fastqc.zip"),
            _path("fastqc_raw", f"{sample_id}_raw_2_fastqc.zip"),
        ])

        if trimming_tool == "trim_galore":
            # Add trim_galore trimming reports
            targets.extend([
                _path("trim", f"{sample_id}_1.fastq.gz_trimming_report.txt"),
                _path("trim", f"{sample_id}_2.fastq.gz_trimming_report.txt"),
            ])
            # Add trim_galore FastQC reports on trimmed reads (ZIP files contain fastqc_data.txt that MultiQC parses)
            targets.extend([
                _path("trim", f"{sample_id}_trimmed_1_fastqc.zip"),
                _path("trim", f"{sample_id}_trimmed_2_fastqc.zip"),
            ])
        else:  # fastp
            # Add FastP quality control JSON
            targets.append(_path("trim", f"{sample_id}.fastp.json"))

    # Conditionally add feature counts summary if enabled
    if is_enabled("featurecounts"):
        targets.append(_path("featurecounts", f"{sample_id}.fc.summary"))

    # Include Salmon metadata so MultiQC waits for quantification.
    if is_enabled("salmon_counts"):
        targets.append(_path("salmon", "aux_info", f"{sample_id}.meta_info.json"))

    # Always add samtools stats files (works for both markdup and regular sorted BAMs)
    targets.extend(
        _path("bam", f"{sample_id}.bam.{ext}")
        for ext in ("stats", "flagstat", "idxstats")
    )

    # Add Picard MarkDuplicates metrics if enabled
    if is_enabled("markduplicates"):
        targets.append(_path("bam", f"{sample_id}.markdup.sorted.MarkDuplicates.metrics.txt"))

    # Add RSeQC outputs if enabled
    targets.extend(get_rseqc_targets(outdir, sample_id, purpose="multiqc"))

    return targets

# Collected target files (fixed typo in function name comment)
def get_target_files(sample_ids):
    """
    Determine target files for the workflow based on enabled modules.
    
    Args:
        sample_ids (list): List of sample IDs to process
        
    Returns:
        list: Paths to all output files required by the workflow
    """
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

        if is_enabled("rseqc"):
            targets.extend(get_rseqc_targets(outdir, sample_id, purpose="workflow"))
            
        if is_enabled("genome_load_keep_memory"):
            targets.append(_path("bam", f"{sample_id}.star_memory_removal.log"))

        # Always include deletion log as final step
        targets.append(_path("logs", f"{sample_id}.deletion.log"))
    
    return targets
