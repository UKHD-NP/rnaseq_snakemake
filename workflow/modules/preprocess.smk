# Rule for Adapter and quality trimming
rule fastp:
    params:
        fastp_params = config['fastp_params']['ffpe']
    input:
        get_paired_fq
    output:
        out1 = os.path.join("{outdir}", "fastp", "{sample_id}_R1.fastp.fastq.gz"),
        out2 = os.path.join("{outdir}", "fastp", "{sample_id}_R2.fastp.fastq.gz"),
        json = os.path.join("{outdir}", "fastp", "{sample_id}.fastp.json"),
        html = os.path.join("{outdir}", "fastp", "{sample_id}.fastp.html"),
        unpaired1 = os.path.join("{outdir}", "fastp", "{sample_id}_R1.fail.fastq.gz"),
        unpaired2 = os.path.join("{outdir}", "fastp", "{sample_id}_R2.fail.fastq.gz") 
    message:
        "{wildcards.sample_id}: Trimming and performing quality control on paired-end FASTQ files"
    log:
        os.path.join("{outdir}", "fastp", "{sample_id}.fastp.log")
    conda:
        "../envs/fastp.yaml"
    threads: 6
    shell: 
        """
        fastp \
            --in1 {input[0]} \
            --in2 {input[1]} \
            --out1 {output.out1} \
            --out2 {output.out2} \
            --json {output.json} \
            --html {output.html} \
            --unpaired1 {output.unpaired1} \
            --unpaired2 {output.unpaired2} \
            --thread {threads} \
            --detect_adapter_for_pe \
            {params.fastp_params} \
            2> {log}
        """

# Rule for Read QC
rule fastqc_before_trimming:
    input:
        get_paired_fq
    output:
        html1 = os.path.join("{outdir}", "fastqc", "{sample_id}_R1.fastqc.html"),
        zip1 = os.path.join("{outdir}", "fastqc", "{sample_id}_R1.fastqc.zip"),
        html2 = os.path.join("{outdir}", "fastqc", "{sample_id}_R2.fastqc.html"),
        zip2 = os.path.join("{outdir}", "fastqc", "{sample_id}_R2.fastqc.zip")
    message:
        "{wildcards.sample_id}: Performing quality control on raw reads"
    log:
        os.path.join("{outdir}", "fastqc", "{sample_id}.fastqc.log")
    conda:
        "../envs/fastqc.yaml"    
    threads: 2
    shell:
        """
        # Run FastQC
        fastqc \
            --quiet \
            --threads {threads} \
            {input[0]} {input[1]} \
            --outdir {wildcards.outdir}/fastqc 2> {log}
        
        fq1_name="$(basename {input[0]})" 
        fq2_name="$(basename {input[1]})"

        # Rename and move output files to match Snakemake expectations
        mv {wildcards.outdir}/fastqc/${{fq1_name%%.*}}_fastqc.html  {output.html1}
        mv {wildcards.outdir}/fastqc/${{fq1_name%%.*}}_fastqc.zip {output.zip1}
        mv {wildcards.outdir}/fastqc/${{fq2_name%%.*}}_fastqc.html  {output.html2}
        mv {wildcards.outdir}/fastqc/${{fq2_name%%.*}}_fastqc.zip {output.zip2}
        """

rule fastqc_after_trimming:
    input:
        get_paired_trimmed_fq
    output:
        html1 = os.path.join("{outdir}", "fastqc", "{sample_id}_R1.fastp.fastqc.html"),
        zip1 = os.path.join("{outdir}", "fastqc", "{sample_id}_R1.fastp.fastqc.zip"),
        html2 = os.path.join("{outdir}", "fastqc", "{sample_id}_R2.fastp.fastqc.html"),
        zip2 = os.path.join("{outdir}", "fastqc", "{sample_id}_R2.fastp.fastqc.zip")
    message:
        "{wildcards.sample_id}: Performing quality control on trimmed reads"
    log:
        os.path.join("{outdir}", "fastqc", "{sample_id}.fastp.fastqc.log")
    conda:
        "../envs/fastqc.yaml"    
    threads: 2
    shell:
        """        
        # Run FastQC
        fastqc \
            --quiet \
            --threads {threads} \
            {input[0]} {input[1]} \
            --outdir {wildcards.outdir}/fastqc 2> {log}

        # Rename and move output files to match Snakemake expectations
        mv {wildcards.outdir}/fastqc/*R1.fastp_fastqc.html  {output.html1}
        mv {wildcards.outdir}/fastqc/*R1.fastp_fastqc.zip {output.zip1}
        mv {wildcards.outdir}/fastqc/*R2.fastp_fastqc.html  {output.html2}
        mv {wildcards.outdir}/fastqc/*R2.fastp_fastqc.zip {output.zip2}
        """
