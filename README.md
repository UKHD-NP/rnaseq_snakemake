# RNA-seq Snakemake Pipeline

Modular RNA-seq workflow built with Snakemake for paired-end short-read data.
The pipeline supports lane merging, optional trimming, optional rRNA removal (SortMeRNA), STAR alignment, duplicate marking, quantification, RSeQC, CPM-normalized bigWig tracks, and per-sample MultiQC reports.

## Table of Contents

- [Workflow Summary](#workflow-summary)
- [Workflow DAG](#workflow-dag)
- [Repository Layout](#repository-layout)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
  - [1. Samplesheet](#1-samplesheet)
  - [2. Reference and modules](#2-reference-and-modules)
- [Local Run](#local-run)
- [Running on HPC with LSF](#running-on-hpc-with-lsf)
  - [Node roles](#node-roles)
  - [Step 1 — Set up Snakemake environment](#step-1---set-up-snakemake-environment)
  - [Step 2 — Clone the pipeline](#step-2---clone-the-pipeline)
  - [Step 3 — Edit configuration](#step-3---edit-configuration)
  - [Step 4 — Update conda-prefix](#step-4---update-conda-prefix-in-the-lsf-profile)
  - [Step 5 — Validate with a dry-run](#step-5---validate-with-a-dry-run)
  - [Step 6 — Submit to HPC](#step-6---submit-to-hpc)
  - [Monitoring jobs](#monitoring-jobs)
- [Output Structure](#output-structure-per-sample)
- [Configuration Reference](#configuration-reference)
  - [Core Settings](#core-settings)
  - [ref Section](#ref-section)
  - [Trimming and Alignment Profiles](#trimming-and-alignment-profiles)
  - [Quantification](#quantification)
  - [Optional Modules](#optional-modules)
  - [RSeQC](#rseqc)
- [MultiQC Notes](#multiqc-notes)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [References](#references)
- [License](#license)

---

## Workflow Summary

For each `sample_id`, the pipeline can run:

1. Merge raw FASTQ lanes (if multiple rows share the same `sample_id`)
2. Trimming (`fastp` or `trim_galore`)
3. FastQC on raw and/or trimmed reads
4. Optional rRNA removal (`SortMeRNA`)
5. STAR genome index generation (if needed) and alignment
6. BAM sorting and optional duplicate marking
7. Quantification (`featureCounts`, `salmon`)
8. Optional modules (`stringtie`, `dupradar`, `arriba`, `RSeQC`, `bigWig`)
9. MultiQC report generation
10. Cleanup of temporary files

## Workflow DAG

Directed Acyclic Graph (DAG) of the pipeline flow (per sample):

```text
samplesheet + reference prep
        |
        v
merge_raw_fastqs
        |
        +--> fastqc_raw
        |
        v
trimming (fastp | trim_galore)
        |
        +--> fastqc_trimmed (trim_galore mode)
        |
        v
sortmerna (optional; removes rRNA reads before alignment)
        |
        v
star_genome_generate (once, if no prebuilt index)
        |
        v
star_align
        |
        v
    sort_bam
        |
        v
mark_duplicates (optional)
        |
        v
    sorted.bam / sorted.markdup.bam
        |
        +--> samtools_stats (flagstat/idxstats/stats)
        |
        +--> featurecounts (optional)
        |
        +--> salmon (optional; uses trimmed FASTQs + tx_fasta)
        |
        +--> stringtie (optional)
        |
        +--> dupradar (optional)
        |
        +--> arriba / fusion (optional; uses FASTA + GTF)
        |
        +--> rseqc (optional; uses BED from gtf2bed)
        |         +--> bam_stat
        |         +--> infer_experiment
        |         +--> inner_distance
        |         +--> read_distribution
        |         +--> read_duplication
        |         +--> read_GC
        |         +--> junction_annotation
        |         +--> junction_saturation
        |         +--> gene_body_coverage
        |         +--> tin
        |
        +--> bigwig (optional; forward + reverse CPM bigWig via deepTools bamCoverage)
        |
        v
     multiqc
        |
        v
    delete_tmp
```

## Repository Layout

```text
rnaseq_snakemake/
  config/config.yml
  workflow/Snakefile
  workflow/modules/*.smk
  workflow/envs/*.yml
  workflow/scripts/*
  test_data/
```

## Requirements

- Linux
- Snakemake ≥ 8 (in a dedicated controller environment)
- Conda/Mamba

> **HPC users:** skip this section and follow [Running on HPC with LSF](#running-on-hpc-with-lsf) instead, which covers environment setup outside your home directory.

For local use, create a minimal controller environment:

```bash
mamba create -n rnaseq_snakemake -c conda-forge -c bioconda snakemake
mamba activate rnaseq_snakemake
```

Each rule uses its own isolated Conda environment defined in `workflow/envs/*.yml`.
Pass `--use-conda` on every Snakemake invocation so these per-rule envs are built and activated automatically.

## Installation

```bash
git clone https://github.com/UKHD-NP/rnaseq_snakemake.git
cd rnaseq_snakemake
```

## Input Files

### 1. Samplesheet

Set `samples_csv` in `config/config.yml` to a CSV with columns:

- `sample_id`
- `fq1`
- `fq2`
- `outdir`

Example:

```csv
sample_id,fq1,fq2,outdir
WT,test_data/raw/SRR6357070_1.fastq.gz,test_data/raw/SRR6357070_2.fastq.gz,test_data/results/WT
WT,test_data/raw/SRR6357071_1.fastq.gz,test_data/raw/SRR6357071_2.fastq.gz,test_data/results/WT
WT,test_data/raw/SRR6357072_1.fastq.gz,test_data/raw/SRR6357072_2.fastq.gz,test_data/results/WT
MUTATION,test_data/raw/SRR6357076_1.fastq.gz,test_data/raw/SRR6357076_2.fastq.gz,test_data/results/MUTATION
```

Notes:
- Repeated `sample_id` rows are treated as lanes and merged before alignment.
- All rows with the same `sample_id` must share the same `outdir`.

### 2. Reference and modules

Edit `config/config.yml`:
- `ref.assembly`: `hg19`, `hg38`, `chm13v2`, `m39`, or `custom`
- Enable/disable optional modules (`featurecounts`, `salmon_counts`, `dupradar`, `fusion`, `rseqc`, etc.)
- Choose a parameter profile (`ffpe`, `total_rna`, etc.)

See [Configuration Reference](#configuration-reference) for all options.

## Local Run

> For cluster execution on HPC, see [Running on HPC with LSF](#running-on-hpc-with-lsf) below.
> The commands here are for single-machine (local) execution only.

**Step 1 — Dry-run first (always).**
Resolves the full DAG and prints every rule that would run — without executing anything:

```bash
snakemake -s workflow/Snakefile --use-conda -n
```

**Step 2 — Optionally verify with the bundled test dataset.**
Runs the full pipeline end-to-end on small test data:

```bash
snakemake -s workflow/Snakefile \
    --configfile config/config_test.yml \
    --use-conda --conda-frontend mamba \
    --cores all
```

**Step 3 — Run with your real config.**

```bash
# Normal run
snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --cores 16

# Rerun only failed/incomplete jobs after fixing an error
snakemake -s workflow/Snakefile --use-conda --conda-frontend mamba --cores 16 --rerun-incomplete
```

> `config/config.yml` is loaded automatically by the Snakefile as the default configfile.
> Pass `--configfile path/to/other.yml` only when you want to override it (e.g. for a test config).

**Re-run MultiQC only for one sample:**

```bash
snakemake -s workflow/Snakefile --use-conda --cores 4 -- results/WT/multiqc/WT.multiqc.html
```

Replace the target path with your sample-specific `outdir`.

## Running on HPC with LSF

This setup uses **IBM Spectrum LSF**.
A ready-made LSF profile is provided at `workflow/profiles/lsf/config.yaml`.

### Node roles

| Node | Purpose | Allowed |
|------|---------|---------|
| `<worker-node>` | Dev, install, testing | ✅ Software install, small runs |
| `<submit-node>` | Job submission only | ✅ Run Snakemake (lightweight), ❌ Processing |
| Cluster nodes | Computation | Jobs submitted automatically via `bsub` |

### Step 1 - Set up Snakemake environment

> **Do this on `<worker-node>`, not on `<submit-node>`.**
> Worker nodes allow software installation. Submission hosts do not.

```bash
ssh YOUR_USERNAME@<worker-node>
```

**Configure conda channels.**
Some HPC clusters ban the `defaults` (Anaconda) channel due to licensing restrictions.
You may need to explicitly restrict to `conda-forge` and `bioconda`:

```bash
cat > ~/.condarc << 'EOF'
channels:
  - conda-forge
  - bioconda
EOF
```

**Load Mamba and initialise your shell.**
This adds `mamba`/`conda` to your `PATH` permanently via `~/.bashrc`:

```bash
module load Mamba   # adjust module name to your site
mamba init bash
source ~/.bashrc   # apply changes to the current shell without re-logging in
```

**Create the Snakemake controller environment outside your home directory.**
Home quota on HPC systems is often limited. Conda environments can easily exceed this — install them on group storage:

```bash
# Set your working directory on group storage (adjust path as needed)
YOUR_WORKDIR="/path/to/group/storage/YOUR_USERNAME"
mkdir -p ${YOUR_WORKDIR}/conda_envs

# Create the controller environment with Snakemake + the LSF executor plugin
mamba create -p ${YOUR_WORKDIR}/conda_envs/snakemake \
    -c conda-forge -c bioconda \
    snakemake \
    snakemake-executor-plugin-lsf \
    -y

# Activate the new environment
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

# Pin numpy/pandas to versions tested with this pipeline's helper scripts
python -m pip install "snakemake==8.*" "snakemake-executor-plugin-lsf" "numpy==1.26.4" "pandas==2.2.3"

# Verify that all three packages are importable and print their versions
python -c "import snakemake, numpy, pandas; print(snakemake.__version__, numpy.__version__, pandas.__version__)"
```

> `snakemake-executor-plugin-lsf` translates Snakemake rule resources (`mem_mb`, `runtime`, `threads`) into `bsub` submission flags automatically — no manual `bsub` scripting needed.

### Step 2 - Clone the pipeline

```bash
cd ${YOUR_WORKDIR}
git clone https://github.com/UKHD-NP/rnaseq_snakemake.git
cd rnaseq_snakemake
```

### Step 3 - Edit configuration

Open `config/config.yml` and set at minimum:
- `samples_csv`: path to your samplesheet CSV
- `ref.assembly`: `hg19`, `hg38`, `chm13v2`, `m39`, or `custom`
- Output directories (via the `outdir` column in the samplesheet)
- Enable/disable optional modules (`featurecounts`, `salmon_counts`, `dupradar`, `fusion`, `rseqc`, etc.)
- Select a parameter profile (`ffpe`, `total_rna`, etc.)

See [Configuration Reference](#configuration-reference) for all options.

### Step 4 - Update `conda-prefix` in the LSF profile

`conda-prefix` tells Snakemake where to build and cache the per-rule conda environments (from `workflow/envs/*.yml`).
All rule environments combined take roughly **5–15 GB** and must live outside your home directory.

Update the placeholder path to your actual working directory:

```bash
sed -i "s|/path/to/group/storage/conda_envs|${YOUR_WORKDIR}/conda_envs|g" \
    workflow/profiles/lsf/config.yaml

# Confirm the replacement was applied correctly
grep "conda-prefix" workflow/profiles/lsf/config.yaml
```

> **Note:** Add the following line to your `~/.bashrc` (once, then `source ~/.bashrc`).
> LSF enforces memory limits per-job, so this variable tells the LSF plugin to
> submit the full `mem_mb` value as a per-job request instead of dividing it per slot:
>
> ```bash
> export SNAKEMAKE_LSF_MEMFMT=perjob
> ```

### Step 5 - Validate with a dry-run

Resolves the full DAG and prints every rule that would run — **without executing or submitting any jobs**.
Always do this before submitting to the cluster to catch config errors, missing inputs, or unexpected rule counts.

```bash
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake
cd ${YOUR_WORKDIR}/rnaseq_snakemake

# Dry-run: prints all rules, checks all inputs, submits nothing
snakemake -s workflow/Snakefile --use-conda -n
```

Confirm that the printed rule count and sample names match expectations before proceeding to Step 6.

> For local testing with the bundled test dataset, see the [Local Run](#local-run) section.

### Step 6 - Submit to HPC

> **Do this on `<submit-node>`**, not on `<worker-node>`.
> Snakemake must run on a submission host to dispatch jobs via `bsub`.

Use `screen` so the Snakemake controller process survives SSH disconnects:

```bash
ssh YOUR_USERNAME@<submit-node>

# Create a named screen session — it keeps running after SSH disconnect
screen -S <session_name>

# Set your working directory (same value as used in Step 1)
YOUR_WORKDIR="/path/to/group/storage/YOUR_USERNAME"

# Activate the Snakemake controller environment
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

# Move into the pipeline directory
cd ${YOUR_WORKDIR}/rnaseq_snakemake

# Launch the pipeline — Snakemake submits each rule as a separate bsub job automatically.
# The config/config.yml is loaded automatically from the Snakefile; no --configfile needed.
# Concurrency is controlled by `jobs:` in workflow/profiles/lsf/config.yaml.
snakemake --profile workflow/profiles/lsf
```

To rerun only failed/incomplete jobs after fixing an error:

```bash
snakemake --profile workflow/profiles/lsf --rerun-incomplete
```

To rerun with the test dataset config:

```bash
snakemake --profile workflow/profiles/lsf --rerun-incomplete --configfile config/config_test.yml
```

Force rerun examples:

```bash
# Force one rule for all matching jobs (e.g. rerun all trim_galore jobs)
snakemake --profile workflow/profiles/lsf --forcerun trim_galore

# Force specific output files (target-level force)
snakemake --profile workflow/profiles/lsf --force \
  test_data/results/SAMPLE_ID/trim/SAMPLE_ID_trimmed_1.fastq.gz \
  test_data/results/SAMPLE_ID/trim/SAMPLE_ID_trimmed_2.fastq.gz

# Force all jobs in the DAG to rerun from scratch
snakemake --profile workflow/profiles/lsf --forceall
```

| `screen` command | Action |
|-----------------|--------|
| `screen -S <session_name>` | Start new named session |
| `Ctrl+A`, then `D` | Detach - session keeps running after SSH disconnect |
| `screen -ls` | List all active sessions |
| `screen -r <session_name>` | Re-attach to session |
| `screen -S <session_name> -X quit` | Kill the named session |

### Monitoring jobs

| `bjobs` command | Action |
|-----------------|--------|
| `bjobs -w` | List all running/pending jobs |
| `bjobs -w -r` | Running only |
| `bjobs -w -p` | Pending only |
| `bjobs -l JOB_ID` | Detailed info for one job |

## Output Structure (Per Sample)

Common outputs in each sample `outdir`:

- `raw_merged/` — merged or symlinked FASTQ files (removed by cleanup when no sample files remain)
- `trim/` — trimmed FASTQ files and trimming reports
- `sortmerna/` — non-rRNA FASTQ files (SortMeRNA output, if enabled)
- `bam/` — STAR-aligned and BAM-derived files
- `featurecounts/` — count matrix and `.fc.summary`
- `salmon/` — quantification outputs
- `rseqc/` — selected RSeQC outputs
- `bigwig/` — CPM-normalized forward/reverse bigWig tracks (if enabled)
- `multiqc/<sample>.multiqc.html`
- `logs/` — rule logs
- `benchmarks/` — Snakemake benchmark files

Final workflow targets are assembled in `rule all` and depend on enabled modules.

## Configuration Reference

Below are the parameters used by the workflow code.

### Core Settings

| Key | Type | Description |
|---|---|---|
| `samples_csv` | string | Path to sample sheet CSV (`sample_id,fq1,fq2,outdir`). |
| `latency-wait` | int | Snakemake filesystem latency wait. |

### `ref` Section

| Key | Type | Description |
|---|---|---|
| `ref.assembly` | string | `hg19`, `hg38`, `chm13v2`, `m39`, or `custom`. |
| `ref.fasta` | string | Required when `assembly: custom`. Can be `.gz`. |
| `ref.gtf` | string | Required when `assembly: custom`. Can be `.gz`. |
| `ref.staridx` | string | Optional prebuilt STAR index directory. |

The runtime keys `ref.tx_fasta` (for quantification) and `ref.bed` (for RSeQC) are generated by `workflow/modules/prepare_genome.smk`.
`ref.bed` is produced via `gtf2bed` in the dedicated env `workflow/envs/gtf2bed.yml` (Perl + gzip/unzip) to keep runs portable across HPC systems.

### Trimming and Alignment Profiles

| Key | Type | Description |
|---|---|---|
| `trimming.enabled` | bool/string/int | Enable trimming-aware branches. |
| `trimming.tool` | string | `fastp` or `trim_galore`. |
| `trimming.param_type` | string | Profile key used in `fastp_params` / `trim_galore_params`. |
| `fastp_params.ffpe` | string | Extra CLI options for fastp FFPE profile. |
| `fastp_params.total_rna` | string | Extra CLI options for fastp total RNA profile. |
| `fastp_params.other` | string | Optional custom profile options. |
| `trim_galore_params.ffpe` | string | Extra CLI options for trim_galore FFPE profile. |
| `trim_galore_params.total_rna` | string | Extra CLI options for trim_galore total RNA profile. |
| `trim_galore_params.other` | string | Optional custom profile options. |
| `alignment.param_type` | string | STAR alignment profile (`ffpe`, `total_rna`, etc.). |
| `star_params.index` | string | STAR genomeGenerate options. |
| `star_params.default` | string | General STAR options profile. |
| `star_params.ffpe` | string | STAR options for FFPE. |
| `star_params.total_rna` | string | STAR options for total RNA. |
| `star_params.fusion` | string | STAR options for Arriba fusion mapping. |
| `genome_load_keep_memory.enabled` | bool/string/int | Enable STAR shared memory cleanup target. **Only set to `true` when `star_params` includes `--genomeLoad LoadAndKeep`.** See note below. |

> **Note on `genome_load_keep_memory`:**
> STAR supports loading the genome index into shared memory (`--genomeLoad LoadAndKeep`) so multiple jobs
> can reuse the same in-memory index instead of reloading it per sample. When this mode is active,
> the genome remains in shared memory after all jobs finish and must be explicitly removed with
> `STAR --genomeLoad Remove`. The `star_remove_shared_memory` rule handles this removal.
>
> Enable `genome_load_keep_memory.enabled: true` **only when** your active `star_params` profile
> includes `--genomeLoad LoadAndKeep`. Otherwise the Remove step will error because no shared genome
> is loaded. In default (`--genomeLoad NoSharedMemory`) mode, keep this setting disabled.
>
> Example `star_params` entry to pair with this setting:
> ```yaml
> star_params:
>   total_rna: "--genomeLoad LoadAndKeep --outSAMtype BAM SortedByCoordinate ..."
> ```
>
> **After the pipeline finishes, verify that shared memory has been released:**
> ```bash
> ipcs -m   # should show no STAR-related shared memory segments
> ```
> If a segment is still listed, remove it manually:
> ```bash
> STAR --genomeLoad Remove --genomeDir /path/to/star/index
> ```

### Quantification

| Key | Type | Description |
|---|---|---|
| `salmon_counts.enabled` | bool/string/int | Enable Salmon quantification outputs. |
| `salmon_counts.param_type` | string | Profile key for `salmon_params`. |
| `salmon_params.ffpe` | string | Salmon CLI options for FFPE profile. |
| `salmon_params.total_rna` | string | Salmon CLI options for total RNA profile. |
| `featurecounts.enabled` | bool/string/int | Enable featureCounts rule and outputs. |
| `featurecounts.feature_type` | string | Optional; default `exon` (`-t` argument). |
| `featurecounts.attribute` | string | Optional; default `gene_id` (`-g` argument). |
| `featurecounts.extra_params` | string | Optional extra featureCounts CLI arguments. |

### Optional Modules

| Key | Type | Description |
|---|---|---|
| `markduplicates.enabled` | bool/string/int | Enable Picard MarkDuplicates branch. |
| `fusion.enabled` | bool/string/int | Enable Arriba fusion calling. |
| `stringtie.enabled` | bool/string/int | Enable StringTie outputs. |
| `dupradar.enabled` | bool/string/int | Enable dupRadar QC. |
| `dupradar.stranded` | int | 0 unstranded, 1 stranded, 2 reverse-stranded. |
| `dupradar.paired` | string | `paired` or `single`. |
| `sortmerna.enabled` | bool/string/int | Enable SortMeRNA rRNA removal (runs after trimming, before alignment/fusion). |
| `sortmerna.delete_sortmerna` | bool/string/int | Delete non-rRNA FASTQ files after pipeline completes (only applies when enabled). |
| `sortmerna.extra_params` | string | Optional extra `sortmerna` CLI arguments. |
| `bigwig.enabled` | bool/string/int | Enable CPM-normalized forward/reverse bigWig generation (deepTools bamCoverage). |
| `bigwig.bin_size` | int | `--binSize` (bp) passed to bamCoverage; smaller = finer resolution, larger file. |
| `bigwig.extra_params` | string | Optional extra `bamCoverage` CLI arguments. |

> **Note on `bigwig`:**
> The rule uses deepTools `--filterRNAstrand forward|reverse`, which assumes a standard reverse-stranded (dUTP) protocol —
> matching this pipeline's default `stringtie.strand: "rf"` / `dupradar.stranded: 2` settings. If your library is
> forward-stranded or unstranded, the `forward`/`reverse` bigWig labels will not match the true transcript strand;
> adjust `workflow/modules/bigwig.smk` accordingly.

### RSeQC

| Key | Type | Description |
|---|---|---|
| `rseqc.enabled` | bool/string/int | Master switch for RSeQC outputs. |
| `rseqc.bam_stat.enabled` | bool | Enable `bam_stat.py`. |
| `rseqc.infer_experiment.enabled` | bool | Enable `infer_experiment.py`. |
| `rseqc.inner_distance.enabled` | bool | Enable `inner_distance.py`. |
| `rseqc.read_distribution.enabled` | bool | Enable `read_distribution.py`. |
| `rseqc.read_duplication.enabled` | bool | Enable `read_duplication.py`. |
| `rseqc.read_GC.enabled` | bool | Enable `read_GC.py`. |
| `rseqc.junction_annotation.enabled` | bool | Enable `junction_annotation.py`. |
| `rseqc.junction_saturation.enabled` | bool | Enable `junction_saturation.py`. |
| `rseqc.gene_body_coverage.enabled` | bool | Enable `geneBody_coverage.py`. |
| `rseqc.tin.enabled` | bool | Enable `tin.py`. |

## MultiQC Notes

- MultiQC config is in `workflow/scripts/multiqc_config.yml`.
- `featurecounts` module ID must stay lowercase in `run_modules` and `sp`.
- MultiQC log for each sample is written to `logs/multiqc/<sample>.multiqc.log`.

## Troubleshooting

- `MissingInputException` in MultiQC:
  - check module toggles and corresponding outputs
  - run dry-run first
- Rule env issues:
  - ensure `--use-conda`
  - delete broken env under `.snakemake/conda/` and rerun
- Large runs:
  - increase `--cores`
  - tune per-rule params in config (aligner/featureCounts/salmon)
- Job killed / out of memory on HPC:
  - check the LSF job log (`bpeek JOB_ID` or `bhist -l JOB_ID`) to confirm out-of-memory (OOM) as the cause
  - **quick fix:** add or increase `mem_mb` for the failing rule in `workflow/profiles/lsf/config.yaml` under `set-resources` — this overrides the rule default without touching the code
  - **permanent fix:** if the rule's default in `workflow/modules/<rule>.smk` under `resources:` is too low, increase `mem_mb` there so the default itself is correct for all runs
- If Snakemake cannot create cache directories in restricted environments, set:

```bash
export XDG_CACHE_HOME=/tmp
```

- If custom references are used, ensure both FASTA and GTF are from the same assembly build.
- For STAR shared-memory issues after a failed run, manually remove the shared memory segment:

```bash
STAR --genomeLoad Remove --genomeDir /path/to/star/index
ipcs -m   # verify no segment remains
```

## Acknowledgments

A huge thank you to Dr. Isabell Bludau, Dr.med.Abigail Suwala, Dr. Paul Kerbs and Quynh Nhu Nguyen from Heidelberg University Hospital and the German Cancer Research Center (DKFZ) for their support, feedback, and contributions to this pipeline.

## References

1. Patel H, Manning J, Ewels P, et al. nf-core/rnaseq [v3.22.2 - Perfect Palladium Penguin]. Zenodo; 2025. https://nf-co.re/rnaseq/3.22.2

## License

Follow the repository [MIT License](MIT_License.md) and tool licenses used in `workflow/envs/`.
