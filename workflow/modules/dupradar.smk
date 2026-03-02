# Rule for Assessment of technical / biological read duplication
rule dupradar:
    input:
        bam = get_bam,
        gtf = config['ref']['gtf']
    output:
        scatter2d = os.path.join("{outdir}", "dupradar", "{sample_id}_duprateExpDens.pdf"),
        boxplot = os.path.join("{outdir}", "dupradar", "{sample_id}_duprateExpBoxplot.pdf"),
        hist = os.path.join("{outdir}", "dupradar", "{sample_id}_expressionHist.pdf"),
        dupmatrix = os.path.join("{outdir}", "dupradar", "{sample_id}_dupMatrix.txt"),
        intercept_slope = os.path.join("{outdir}", "dupradar", "{sample_id}_intercept_slope.txt")
    params:
        prefix = lambda wildcards: os.path.join(wildcards.outdir, "dupradar", wildcards.sample_id),
        stranded = config.get('dupradar', {}).get('stranded', 2),  # 0=unstranded, 1=stranded, 2=reverse
        paired = config.get('dupradar', {}).get('paired', 'paired'),  # 'paired' or 'single'
        script = os.path.join(workflow.basedir, "scripts", "dupradar.r")
    conda:
        os.path.join(workflow.basedir, "envs", "dupradar.yml")
    message:
        "{wildcards.sample_id}: Running dupRadar to evaluate technical and biological read duplication"
    threads: 4
    resources:
        mem_mb = 16384
    log:
        os.path.join("{outdir}", "logs", "dupradar", "{sample_id}.dupradar.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "dupradar.{sample_id}.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.scatter2d})
        mkdir -p $(dirname {log})

        if [ -f "{params.script}" ]; then
            Rscript "{params.script}" \
                {input.bam} \
                {params.prefix} \
                {input.gtf} \
                {params.stranded} \
                {params.paired} \
                {threads} 2> {log} || {{ echo "[ERROR] dupRadar failed." >> {log}; exit 1; }}
        else
            echo "[ERROR] dupradar.r script not found: {params.script}" > {log}
            exit 1
        fi
        """
