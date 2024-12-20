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
            os.path.join(trim_dir, wildcards.sample_id + "_R1.fastp.fastq.gz"),
            os.path.join(trim_dir, wildcards.sample_id + "_R2.fastp.fastq.gz")
            ]
    else:
        sample_data = samplesheet[samplesheet['sample_id'] == wildcards.sample_id].iloc[0]
        return [
            sample_data['fq1'],
            sample_data['fq2']
            ]


def get_input_multiqc(wildcards):
    outdir = get_outdir(wildcards.sample_id)      
    targets = []

    # Fastp output
    targets += [
        os.path.join(outdir, "fastp" , f"{wildcards.sample_id}.fastp.json")
    ]

    # STAR output
    targets += [
        os.path.join(outdir, "bam", f"{wildcards.sample_id}.Log.final.out")
    ]

    # featureCounts output
    targets += [
        os.path.join(outdir, "featurecounts", f"{wildcards.sample_id}.fc.summary")
    ]
    
    return targets


def get_target_files(sample_ids):
    multiqc_files = [os.path.join(get_outdir(sample_id), "multiqc", f"{sample_id}.multiqc.html") for sample_id in sample_ids]
    salmon_files = [os.path.join(get_outdir(sample_id), "salmon") for sample_id in sample_ids]
    fc_files = [os.path.join(get_outdir(sample_id), "featurecounts", f"{sample_id}.fc") for sample_id in sample_ids]
    fusion_files = [os.path.join(get_outdir(sample_id), "arriba", f"{sample_id}.fusions.tsv") for sample_id in sample_ids]
    
    targets = multiqc_files + salmon_files + fc_files + fusion_files
    return targets
