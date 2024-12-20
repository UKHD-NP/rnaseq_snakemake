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
        "../envs/fastp.yml"
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
            --report_title "{wildcards.sample_id} fastp report" \
            {params.fastp_params} \
            2>{log}
        """

rule multiqc:
    params:
        multiqc_cfg = config['multiqc_cfg'],
        outdir = lambda w: os.path.join(f"{w.outdir}", "multiqc")
    input:
        get_input_multiqc
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
            --config {params.multiqc_cfg} \
            2>{log}
        """
