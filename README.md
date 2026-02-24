# RNA-seq Snakemake Pipeline

Modular RNA-seq workflow built with Snakemake for paired-end short-read data.  
The pipeline supports lane merging, optional trimming, STAR alignment, duplicate marking, quantification, RSeQC, and per-sample MultiQC reports.

## Workflow Summary

For each `sample_id`, the pipeline can run:

1. Merge raw FASTQ lanes (if multiple rows share the same `sample_id`)
2. Trimming (`fastp` or `trim_galore`)
3. FastQC on raw and/or trimmed reads
4. STAR genome index generation (if needed) and alignment
5. BAM sorting and optional duplicate marking
6. Quantification (`featureCounts`, `salmon`)
7. Optional modules (`stringtie`, `dupradar`, `arriba`, `RSeQC`)
8. MultiQC report generation
9. Cleanup of temporary files

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

- Linux environment
- Snakemake (recommended: recent 8.x/9.x)
- Conda or Mamba/Micromamba for per-rule environments (`--use-conda`)

## Quick Start

### 1. Prepare the sample sheet

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

- Repeated `sample_id` is allowed for lane merging.
- All rows with the same `sample_id` must share the same `outdir`.

### 2. Configure references and modules

Edit `config/config.yml`:

- select `ref.assembly` (`hg19`, `hg38`, `chm13v2`, `m39`, `custom`)
- enable/disable modules as needed
- choose parameter profiles (`ffpe`, `total_rna`, etc.)

### 3. Dry-run

```bash
snakemake -n -s workflow/Snakefile --configfile config/config.yml
```

### 4. Run

```bash
snakemake --use-conda --cores 16 -s workflow/Snakefile --configfile config/config.yml
```

### 5. Re-run only MultiQC for one sample

```bash
snakemake --use-conda --cores 4 -s workflow/Snakefile --configfile config/config.yml -- results/WT/multiqc/WT.multiqc.html
```

Replace the target path with your sample-specific `outdir`.

## Running on DKFZ HPC (LSF)

The DKFZ cluster uses **IBM Spectrum LSF** (not SLURM).
A ready-made LSF profile is provided at `workflow/profiles/lsf/config.yaml`.

### Node roles at DKFZ

| Node | Purpose | Allowed |
|------|---------|---------|
| `odcf-worker01/02` | Dev, install, testing | ✅ Software install, small runs |
| `bsub01/02` | Job submission only | ✅ Run Snakemake (lightweight), ❌ Processing |
| Cluster nodes | Computation | Jobs submitted automatically via `bsub` |

### Step 1 — Set up Snakemake environment (on odcf-worker01)

Worker nodes allow software installation; submission hosts (`bsub01`) do not.

```bash
ssh YOUR_USERNAME@odcf-worker01.dkfz.de
```

Configure conda channels — required by DKFZ because the `defaults` channel (Anaconda) is banned due to licensing:

```bash
cat > ~/.condarc << 'EOF'
channels:
  - conda-forge
  - bioconda
EOF
```

Load Mamba and initialise your shell:

```bash
module load Mamba/24.11.2-1
mamba init bash
source ~/.bashrc
```

Create the Snakemake controller environment **outside home** (home quota is only 20 GB):

```bash
YOUR_WORKDIR="/omics/groups/OE0146/internal/YOUR_USERNAME"
mkdir -p ${YOUR_WORKDIR}/conda_envs

mamba create -p ${YOUR_WORKDIR}/conda_envs/snakemake \
    -c conda-forge -c bioconda \
    snakemake \
    snakemake-executor-plugin-lsf \
    -y
```

`snakemake-executor-plugin-lsf` lets Snakemake translate rule resources (`mem_mb`, `runtime`, `threads`) into `bsub` flags automatically.

### Step 2 — Clone the pipeline

```bash
cd ${YOUR_WORKDIR}
git clone https://github.com/UKHD-NPS/rnaseq_snakemake.git
cd rnaseq_snakemake
```

### Step 3 — Edit configuration

Edit `config/config.yml`: set `samples_csv`, `ref.assembly`, output directories, and enable/disable modules.

### Step 4 — Update conda-prefix in the LSF profile

Open `workflow/profiles/lsf/config.yaml` and update the `conda-prefix` line, or use sed:

```bash
sed -i "s|/omics/odcf/analysis/YOUR_GROUP/conda_envs|${YOUR_WORKDIR}/conda_envs|g" \
    workflow/profiles/lsf/config.yaml
```

**Why this matters:** `conda-prefix` tells Snakemake where to build and cache per-rule conda environments (from `workflow/envs/*.yml`). All rule environments together take 5–15 GB. This path must be outside home to avoid hitting the 20 GB home quota.

Verify it was applied:

```bash
grep "conda-prefix" workflow/profiles/lsf/config.yaml
```

### Step 5 — Dry-run (validate without submitting any jobs)

```bash
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

snakemake -s workflow/Snakefile \
    --configfile config/config.yml \
    --use-conda -n
```

A dry-run prints every rule Snakemake would execute without running anything. Confirm the job count and sample names look correct before submitting to the cluster.

### Step 6 — Run on HPC (from bsub01)

Snakemake must be launched from a **submission host** (`bsub01` or `bsub02`).
Use `screen` to keep the session alive if SSH disconnects.

```bash
ssh YOUR_USERNAME@bsub01.lsf.dkfz.de

# Activate Snakemake env
module load Mamba/24.11.2-1
mamba activate ${YOUR_WORKDIR}/conda_envs/snakemake

# Go to the pipeline directory
cd ${YOUR_WORKDIR}/rnaseq_snakemake

# Start a persistent screen session (survives SSH disconnect)
screen -S rnaseq

# Run — Snakemake submits each rule as a bsub job automatically
snakemake --profile workflow/profiles/lsf -j 100
```

| `screen` command | Action |
|-----------------|--------|
| `screen -S rnaseq` | Start new named session |
| `Ctrl+A`, then `D` | Detach — session keeps running after SSH disconnect |
| `screen -ls` | List all active sessions |
| `screen -r rnaseq` | Re-attach to session |

### Monitoring jobs

```bash
bjobs -w           # list all running/pending jobs
bjobs -w -r        # running only
bjobs -w -p        # pending only
bjobs -l JOB_ID    # detailed info for one job
```

## Output Structure (Per Sample)

Common outputs in each sample `outdir`:

- `raw_merged/` merged or symlinked FASTQ files (removed by cleanup when no sample files remain)
- `trim/` trimmed FASTQ files and trimming reports
- `bam/` STAR and BAM-derived files
- `featurecounts/` count matrix and `.fc.summary`
- `salmon/` quantification outputs
- `rseqc/` selected RSeQC outputs
- `multiqc/<sample>.multiqc.html`
- `logs/` rule logs
- `benchmarks/` Snakemake benchmark files

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
| `genome_load_keep_memory.enabled` | bool/string/int | Enable STAR shared memory cleanup target. |

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
| `featurecounts.params` | string | Optional extra featureCounts CLI arguments. |

### Optional Modules

| Key | Type | Description |
|---|---|---|
| `markduplicates.enabled` | bool/string/int | Enable Picard MarkDuplicates branch. |
| `fusion.enabled` | bool/string/int | Enable Arriba fusion calling. |
| `stringtie.enabled` | bool/string/int | Enable StringTie outputs. |
| `dupradar.enabled` | bool/string/int | Enable dupRadar QC. |
| `dupradar.stranded` | int | 0 unstranded, 1 stranded, 2 reverse-stranded. |
| `dupradar.paired` | string | `paired` or `single`. |

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

- MultiQC config is in `workflow/scripts/multiqc_config.yaml`.
- `featurecounts` module ID must stay lowercase in `run_modules` and `sp`.
- MultiQC log for each sample is written to `logs/multiqc/<sample>.multiqc.log`.

## Troubleshooting

- If Snakemake cannot create cache directories in restricted environments, set:

```bash
export XDG_CACHE_HOME=/tmp
```

- If custom references are used, ensure both FASTA and GTF are from the same assembly build.

## Acknowledgments

A huge thank you to Dr. Isabell Bludau, Dr.med.Abigail Suwala, Dr. Paul Kerbs and Quynh Nhu Nguyen from Heidelberg University Hospital and the German Cancer Research Center (DKFZ) for their support, feedback, and contributions to this pipeline.

## License

Follow the repository [MIT License](MIT_License.md) and tool licenses used in `workflow/envs/`.
