def is_paired_end(sample_id):
    # Select the relevant row(s) for the given sample_name
    sample_units = samplesheet[samplesheet['sample_id'] == sample_id]
    
    # Check if 'fq2' is null (i.e., missing) for all the rows corresponding to this sample
    fq2_null = sample_units["fq2"].isnull()
    
    # Determine if this sample is paired-end
    paired = ~fq2_null
    
    # Check if all associated rows are either paired-end or single-end
    all_paired = paired.all()
    all_single = (~paired).all()
    
    # Assert that the sample is either completely paired-end or single-end
    assert (
        all_single or all_paired
    ), f"Invalid units for sample {sample_id}, must be all paired-end or all single-end."
    
    return all_paired


def get_outdir(sample_id):
    sample_data = samplesheet[samplesheet['sample_id'] == sample_id].iloc[0]
    return sample_data['outdir']


def get_paired_fq(wildcards):
    if is_paired_end(wildcards.sample_id):
        sample_data = samplesheet[samplesheet['sample_id'] == wildcards.sample_id].iloc[0]
        return [sample_data['fq1'], sample_data['fq2']]  # Return a list of file paths
    else:
        raise ValueError(f"Sample {wildcards.sample} is not paired-end.")


def get_paired_trimmed_fq(wildcards):
    if is_paired_end(wildcards.sample_id):
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
    else:
        raise ValueError(f"Sample {wildcards.sample_id} is not paired-end.")


