# Get the output directory for a given sample id
def get_outdir(sample_id):
    sample_data = samplesheet[samplesheet['sample_id'] == sample_id].iloc[0]
    return sample_data['outdir']


# Get list of paths to raw fastq files from samplesheet
def get_paired_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_id'] == wildcards.sample_id].iloc[0]
    return [
        sample_data['fq1'],
        sample_data['fq2']
        ]


# Get list of paths to trimmed fastq files (if trimming is disabled, return raw fastq files)
def get_paired_trimmed_fq(wildcards):
    if config['trimming']['enabled']:
        trim_dir = os.path.join(get_outdir(wildcards.sample_id), "fastp")
        return [
            os.path.join(trim_dir, wildcards.sample_id + "_1.fastp.fastq.gz"),
            os.path.join(trim_dir, wildcards.sample_id + "_2.fastp.fastq.gz")
            ]
    else:
        sample_data = samplesheet[samplesheet['sample_id'] == wildcards.sample_id].iloc[0]
        return [
            sample_data['fq1'],
            sample_data['fq2']
            ]


def get_input_multiqc_sample(wildcards):
    outdir = get_outdir(wildcards.sample_id)      
    
    # FastQC output
    targets = [
        os.path.join(get_outdir(wildcards.sample_id), "fastqc" , f"{wildcards.sample_id}_1.fastqc.html"),
        os.path.join(get_outdir(wildcards.sample_id), "fastqc" , f"{wildcards.sample_id}_2.fastqc.html"),
        os.path.join(get_outdir(wildcards.sample_id), "fastqc" , f"{wildcards.sample_id}_1.fastp.fastqc.html"),
        os.path.join(get_outdir(wildcards.sample_id), "fastqc" , f"{wildcards.sample_id}_2.fastp.fastqc.html")
        ]
    
    # Fastp output
    targets += [
        os.path.join(get_outdir(wildcards.sample_id), "fastp" , f"{wildcards.sample_id}.fastp.json")
        ]
    
    # Samtools stats
    targets += [
        os.path.join(get_outdir(wildcards.sample_id), "samtools_stats" , f"{wildcards.sample_id}.samtools.stats"),
        os.path.join(get_outdir(wildcards.sample_id), "samtools_stats" , f"{wildcards.sample_id}.samtools.flagstats"),
        os.path.join(get_outdir(wildcards.sample_id), "samtools_stats" , f"{wildcards.sample_id}.samtools.idxstats")
        ]

    return targets


def get_target_files(sample_ids):
    multiqc_files = [os.path.join(get_outdir(sample_id), "multiqc", f"{sample_id}.multiqc.html") for sample_id in sample_ids]
    targets = multiqc_files
    return targets
