rule salmon:
    input:
        tx_bam = os.path.join("{outdir}", "{sample_id}", "bam", "{sample_id}.tx.bam"),
        tx_fasta = config['ref']['tx_fasta'],
        gtf = config['ref']['gtf']
    output:
        directory(os.path.join("{outdir}", "{sample_id}", "salmon"))
    message:
        "{wildcards.sample_id}: Count transcripts with Salmon"
    log:
        os.path.join("{outdir}", "{sample_id}", "salmon", "{sample_id}.salmon.log")
    conda:
        "../envs/salmon.yml"
    threads: 12
    shell:
        """
        salmon quant \
            --threads {threads}
            --targets {input.tx_fasta} \
            --alignments {input.tx_bam} \
            --libType A \
            --geneMap {input.gtf} \
            --minAssignedFrags 5 \
            --output {output} 2> {log}
        """
