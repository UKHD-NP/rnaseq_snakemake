rule bamcoverage_bigwig:
    # Generate strand-specific CPM-normalized bigWig tracks for IGV (reverse-stranded/dUTP protocol)
    input:
        bam = get_bam,
        bai = get_bam_bai
    output:
        bigwig = os.path.join("{outdir}", "bigwig", "{sample_id}.{strand}.CPM.bw")
    wildcard_constraints:
        strand = "forward|reverse"
    params:
        bin_size = config.get('bigwig', {}).get('bin_size', 10),
        extra_params = config.get('bigwig', {}).get('extra_params', '')
    conda:
        os.path.join(workflow.basedir, "envs", "deeptools.yml")
    message:
        "{wildcards.sample_id}: Generating {wildcards.strand}-strand CPM bigWig"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 4096 + (attempt - 1) * 2048,
        runtime = lambda wildcards, attempt: attempt * 480
    log:
        os.path.join("{outdir}", "logs", "bigwig", "{sample_id}.{strand}.bamCoverage.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.{strand}.bamCoverage.benchmark.txt")
    shell:
        """
        ulimit -Sn $(ulimit -Hn) 2>/dev/null || ulimit -n 65536 2>/dev/null || true
        echo "[INFO] File descriptor limit: $(ulimit -n)" >> "{log}"

        mkdir -p $(dirname {output.bigwig})
        mkdir -p $(dirname {log})

        echo "[INFO] bamCoverage start: $(date)" >> "{log}"

        bamCoverage \
            --numberOfProcessors {threads} \
            --binSize {params.bin_size} \
            --normalizeUsing CPM \
            --bam {input.bam} \
            --filterRNAstrand {wildcards.strand} \
            --outFileName {output.bigwig} \
            --outFileFormat bigwig \
            {params.extra_params} \
            > {log} 2>&1 || {{ echo "[ERROR] bamCoverage ({wildcards.strand}) failed." >> {log}; exit 1; }}

        if [ ! -s {output.bigwig} ]; then
            echo "[ERROR] bigWig output is missing or empty: {output.bigwig}" >> {log}
            exit 1
        fi
        """
