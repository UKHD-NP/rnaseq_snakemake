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
            --threads {threads} \
            --targets {input.tx_fasta} \
            --alignments {input.tx_bam} \
            --gcBias --seqBias \
            --libType A \
            --geneMap {input.gtf} \
            --minAssignedFrags 5 \
            --output {output} 2>{log}
        
        mv {output}/quant.sf {output}/{wildcards.sample_id}.quant.sf
        mv {output}/quant.genes.sf {output}/{wildcards.sample_id}.quant.genes.sf
        """

rule featurecounts:
    input:
        bam = os.path.join("{outdir}", "{sample_id}", "bam", "{sample_id}.bam"),
        gtf = config['ref']['gtf'],
        fasta = config['ref']['fasta']
    output:
        os.path.join("{outdir}", "{sample_id}", "featurecounts","{sample_id}.fc.txt")
    message:
        "{wildcards.sample_id}: Count reads with featureCounts (Subread)"
    log:
        os.path.join("{outdir}", "{sample_id}", "featurecounts", "{sample_id}.fc.log")
    conda:
        "../envs/subread.yml"
    threads: 8
    shell:
        """
        featureCounts \
            -a {input.gtf} \
            -o {output} \
            -T {threads} \
            -p \
            --countReadPairs \
            -t exon \
            -g gene_id \
            -Q 20 \
            -J \
            -G {input.fasta} \
            {input.bam} >{log} 2>&1
        """
