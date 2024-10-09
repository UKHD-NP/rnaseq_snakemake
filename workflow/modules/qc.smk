rule samtools_stats:
    input:
        bam = os.path.join("{outdir}", "bam" , "{sample_id}.bam"),
        ref_fasta = config['ref']['fasta']
    output:
        stats = os.path.join("{outdir}", "samtools_stats" , "{sample_id}.samtools.stats"),
        flagstat = os.path.join("{outdir}", "samtools_stats" , "{sample_id}.samtools.flagstats"),
        idxstats = os.path.join("{outdir}", "samtools_stats" , "{sample_id}.samtools.idxstats")   
    message:
        "{wildcards.sample_id}: Generates SAMtools statistics."
    log:
        os.path.join("{outdir}", "samtools_stats" , "{sample_id}.samtools.log")
    conda:
        "../envs/samtools.yml"
    threads: 3
    shell:
        """
         # Generates SAMtools statistics.
        samtools stats \
            --threads 1 \
            --reference {input.ref_fasta} \
            {input.bam} \
            > {output.stats} 2>> {log} &
        
        # Counts the number of alignments in the BAM file for each FLAG type
        samtools flagstat \
            --threads 1 \
            {input.bam} \
            > {output.flagstat} 2>> {log} &
        
        # Produces index statistics for the sorted BAM file
        samtools idxstats \
            {input.bam} \
            > {output.idxstats} 2>> {log}
        """

rule multiqc_per_sample:
    params:
        multiqc_cfg = config['multiqc_cfg'],
        outdir = lambda w: os.path.join(f"{w.outdir}", "multiqc")
    input:
        get_input_multiqc_sample
    output:
        os.path.join("{outdir}", "multiqc", "{sample_id}.multiqc.html")
    log:
        os.path.join("{outdir}", "multiqc", "{sample_id}.multiqc.log")
    conda:
        "../envs/multiqc.yml"
    threads: 1
    message: 
        "{wildcards.sample_id}: Running MultiQC."
    shell:
        """
        multiqc {input} \
            --outdir {params.outdir} \
            --filename {wildcards.sample_id}.multiqc.html \
            --force \
            --config {params.multiqc_cfg} 2> {log}
        """

rule multiqc_all_samples:
    params:
        multiqc_cfg = config['multiqc_cfg'],
        outdir = lambda w: os.path.join(f"{w.outdir}", "multiqc")
    input:
        get_input_multiqc_all
    output:
        os.path.join("{outdir}", "multiqc", "all.multiqc.html")
    log:
        os.path.join("{outdir}", "multiqc", "all.multiqc.log")
    conda:
        "../envs/multiqc.yml"
    threads: 1
    message: 
        "All samples: Running MultiQC."
    shell:
        """
        multiqc {input} \
            --outdir {params.outdir} \
            --filename {wildcards.sample_id}.multiqc.html \
            --force \
            --config {params.multiqc_cfg} 2> {log}
        """
