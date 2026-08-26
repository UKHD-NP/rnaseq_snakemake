rule ribodetector:
    # Remove ribosomal RNA (rRNA) reads with RiboDetector (deep-learning based, no reference DB needed)
    input:
        fq = get_paired_trimmed_fq
    output:
        fq1 = os.path.join("{outdir}", "ribodetector", "{sample_id}.nonrna_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "ribodetector", "{sample_id}.nonrna_2.fastq.gz"),
        summary = os.path.join("{outdir}", "ribodetector", "{sample_id}.ribodetector.log")
    params:
        read_length = config.get('ribodetector', {}).get('read_length', 100),
        chunk_size = config.get('ribodetector', {}).get('chunk_size', 256),
        extra_params = config.get('ribodetector', {}).get('extra_params', '')
    conda:
        os.path.join(workflow.basedir, "envs", "ribodetector.yml")
    message:
        "{wildcards.sample_id}: Removing rRNA reads with RiboDetector"
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32768,
        runtime = lambda wildcards, attempt: attempt * 240
    log:
        os.path.join("{outdir}", "logs", "ribodetector", "{sample_id}.ribodetector.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "ribodetector.{sample_id}.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.fq1})
        mkdir -p $(dirname {log})

        ribodetector_cpu \
            -i {input.fq[0]} {input.fq[1]} \
            -o {output.fq1} {output.fq2} \
            -l {params.read_length} \
            -t {threads} \
            --chunk_size {params.chunk_size} \
            --log {output.summary} \
            {params.extra_params} \
            > {log} 2>&1 || {{ echo "[ERROR] RiboDetector failed." >> {log}; exit 1; }}

        if [ ! -s {output.fq1} ] || [ ! -s {output.fq2} ]; then
            echo "[ERROR] RiboDetector non-rRNA FASTQ outputs are missing or empty." >> {log}
            exit 1
        fi
        """
