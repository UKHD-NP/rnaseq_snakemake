rule delete_fastqs:
    input:
        fq1 = os.path.join("{outdir}", "fastp", "{sample_id}_R1.fastp.fastq.gz"),
        fq2 = os.path.join("{outdir}", "fastp", "{sample_id}_R2.fastp.fastq.gz"),
        fq1_fail = os.path.join("{outdir}", "fastp", "{sample_id}_R1.fail.fastq.gz"),
        fq2_fail = os.path.join("{outdir}", "fastp", "{sample_id}_R2.fail.fastq.gz"),
        fc_sum = os.path.join("{outdir}", "featurecounts","{sample_id}.fc.summary"),
        ar_tsv = os.path.join("{outdir}", "arriba", "{sample_id}.fusions.tsv")
    output:
        os.path.join("{outdir}", "fastp", "{sample_id}.fastqdeletion.log")
    threads: 1
    shell:
        """
        rm {input.fq1} {input.fq2} {input.fq1_fail} {input.fq2_fail} > {output}
        """
