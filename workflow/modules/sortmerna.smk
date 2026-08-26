SORTMERNA_IDX_DIR = os.path.join(sortmerna_db_dir, "idx")

rule sortmerna_index:
    # Build the shared SortMeRNA index from the rRNA reference database (built once, reused by all samples)
    input:
        config["sortmerna"]["ref_db"]
    output:
        touch(os.path.join(SORTMERNA_IDX_DIR, ".sortmerna_idx.done"))
    params:
        refs = " ".join(f"--ref {f}" for f in config["sortmerna"]["ref_db"]),
        idx_dir = SORTMERNA_IDX_DIR,
        workdir = os.path.join(sortmerna_db_dir, "idx_workdir")
    conda:
        os.path.join(workflow.basedir, "envs", "sortmerna.yml")
    message:
        "Building SortMeRNA index from rRNA reference database"
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8192,
        runtime = lambda wildcards, attempt: attempt * 120
    log:
        os.path.join("references", "logs", "sortmerna_index.log")
    shell:
        """
        mkdir -p {params.idx_dir}
        mkdir -p $(dirname {log})
        rm -rf {params.workdir}

        sortmerna \
            {params.refs} \
            --idx-dir {params.idx_dir} \
            --workdir {params.workdir} \
            --task 5 \
            --threads {threads} \
            -v \
            > {log} 2>&1 || {{ echo "[ERROR] SortMeRNA indexing failed." >> {log}; exit 1; }}

        rm -rf {params.workdir}
        """

rule sortmerna:
    # Remove ribosomal RNA (rRNA) reads with SortMeRNA (reference-based, CPU-only)
    input:
        fq = get_paired_trimmed_fq,
        idx_done = os.path.join(SORTMERNA_IDX_DIR, ".sortmerna_idx.done")
    output:
        fq1 = os.path.join("{outdir}", "sortmerna", "{sample_id}.nonrna_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "sortmerna", "{sample_id}.nonrna_2.fastq.gz"),
        summary = os.path.join("{outdir}", "sortmerna", "{sample_id}.sortmerna.log")
    params:
        refs = " ".join(f"--ref {f}" for f in config["sortmerna"]["ref_db"]),
        idx_dir = SORTMERNA_IDX_DIR,
        workdir = os.path.join("{outdir}", "sortmerna", "wd_{sample_id}"),
        other_prefix = os.path.join("{outdir}", "sortmerna", "{sample_id}.nonrna"),
        aligned_prefix = os.path.join("{outdir}", "sortmerna", "{sample_id}.rrna"),
        extra_params = config.get('sortmerna', {}).get('extra_params', '')
    conda:
        os.path.join(workflow.basedir, "envs", "sortmerna.yml")
    message:
        "{wildcards.sample_id}: Removing rRNA reads with SortMeRNA"
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16384,
        runtime = lambda wildcards, attempt: attempt * 480
    log:
        os.path.join("{outdir}", "logs", "sortmerna", "{sample_id}.sortmerna.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "sortmerna.{sample_id}.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.fq1})
        mkdir -p $(dirname {log})
        rm -rf {params.workdir}
        mkdir -p {params.workdir}

        sortmerna \
            {params.refs} \
            --idx-dir {params.idx_dir} \
            --workdir {params.workdir} \
            --reads {input.fq[0]} \
            --reads {input.fq[1]} \
            --paired_in \
            --out2 \
            --fastx \
            --other {params.other_prefix} \
            --aligned {params.aligned_prefix} \
            --zip-out yes \
            --threads {threads} \
            -v \
            {params.extra_params} \
            > {log} 2>&1 || {{ echo "[ERROR] SortMeRNA failed." >> {log}; exit 1; }}

        mv {params.other_prefix}_fwd.fq.gz {output.fq1} 2>>{log} || {{ echo "[ERROR] Missing SortMeRNA non-rRNA R1 output." >> {log}; exit 1; }}
        mv {params.other_prefix}_rev.fq.gz {output.fq2} 2>>{log} || {{ echo "[ERROR] Missing SortMeRNA non-rRNA R2 output." >> {log}; exit 1; }}
        mv {params.aligned_prefix}.log {output.summary} 2>>{log} || {{ echo "[ERROR] Missing SortMeRNA summary log." >> {log}; exit 1; }}

        rm -rf {params.workdir} {params.aligned_prefix}_fwd.fq.gz {params.aligned_prefix}_rev.fq.gz

        if [ ! -s {output.fq1} ] || [ ! -s {output.fq2} ]; then
            echo "[ERROR] SortMeRNA non-rRNA FASTQ outputs are missing or empty." >> {log}
            exit 1
        fi
        """
