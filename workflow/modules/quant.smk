rule salmon:
    # Count transcripts with Salmon
    input:
        tx_bam = os.path.join("{outdir}", "bam", "{sample_id}.tx.bam"),
        tx_fasta = config['ref']['tx_fasta'],
        gtf = config['ref']['gtf']
    output:
        quant_sf = os.path.join("{outdir}", "salmon", "{sample_id}.quant.sf"),
        quant_genes_sf = os.path.join("{outdir}", "salmon", "{sample_id}.quant.genes.sf"),
        meta_info = os.path.join("{outdir}", "salmon", "aux_info", "{sample_id}.meta_info.json")
    params:
        other_params = lambda wildcards: config["salmon_params"].get(
            config.get("salmon_counts", {}).get("param_type", "total_rna"),
            ""
        ),
        outdir = lambda wildcards: os.path.join(wildcards.outdir, "salmon"),
    conda:
        os.path.join(workflow.basedir, "envs", "salmon.yml")
    message:
        "{wildcards.sample_id}: Count transcripts with Salmon"
    threads: 12
    resources:
        mem_mb = 36864
    log:
        os.path.join("{outdir}", "logs", "salmon", "{sample_id}.salmon.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "salmon.{sample_id}.benchmark.txt")
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        salmon quant \
            --threads {threads} \
            --targets {input.tx_fasta} \
            --alignments {input.tx_bam} \
            --geneMap {input.gtf} \
            --minAssignedFrags 5 \
            --output {params.outdir} \
            {params.other_params} 2> {log} || {{ echo "[ERROR] Salmon quantification failed." >> {log}; exit 1; }}

        # Rename output files for clarity
        mv {params.outdir}/quant.sf {output.quant_sf} || {{ echo "[ERROR] Missing Salmon output: quant.sf." >> {log}; exit 1; }}
        mv {params.outdir}/quant.genes.sf {output.quant_genes_sf} || {{ echo "[ERROR] Missing Salmon output: quant.genes.sf." >> {log}; exit 1; }}
        cp {params.outdir}/aux_info/meta_info.json {output.meta_info} || {{ echo "[ERROR] Missing Salmon output: aux_info/meta_info.json." >> {log}; exit 1; }}
        """

rule featurecounts:
    # Count reads with featureCounts (Subread)
    input:
        bam = get_bam,
        gtf = config['ref']['gtf'],
        fasta = config['ref']['fasta']
    output:
        counts = os.path.join("{outdir}", "featurecounts", "{sample_id}.fc"),
        summary = os.path.join("{outdir}", "featurecounts", "{sample_id}.fc.summary"),
        jcounts = os.path.join("{outdir}", "featurecounts", "{sample_id}.fc.jcounts")
    params:
        # Additional parameters from config
        extra_params = config.get('featurecounts', {}).get('extra_params', ""),
        # Define feature type and attribute from config with defaults
        feature_type = config.get('featurecounts', {}).get('feature_type', "exon"),
        attribute = config.get('featurecounts', {}).get('attribute', "gene_id"),
    conda:
        os.path.join(workflow.basedir, "envs", "subread.yml")
    message:
        "{wildcards.sample_id}: Count reads with featureCounts (Subread)"
    threads: 6
    resources:
        mem_mb = 36864
    log:
        os.path.join("{outdir}", "logs", "featurecounts", "{sample_id}.fc.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "featurecounts.{sample_id}.benchmark.txt")
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.counts})
        mkdir -p $(dirname {log})

        # Log which BAM file is being used
        echo "[INFO] Using BAM file: {input.bam}" > {log}
        echo "[INFO] Feature type: {params.feature_type}" >> {log}
        echo "[INFO] Attribute: {params.attribute}" >> {log}

        # Export memory limit for better resource usage
        export MALLOC_ARENA_MAX=4

        featureCounts \
            -a {input.gtf} \
            -o {output.counts} \
            -T {threads} \
            --countReadPairs \
            -p \
            -t {params.feature_type} \
            -g {params.attribute} \
            -Q 20 \
            -J \
            -G {input.fasta} \
            {params.extra_params} \
            {input.bam} >> {log} 2>&1 || {{ echo "[ERROR] featureCounts failed." >> {log}; exit 1; }}
        """
