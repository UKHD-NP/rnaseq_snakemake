config.setdefault("ref", {})

STAR_INDEX_DIR = str(config["ref"].get("star_index", "")).strip()
if not STAR_INDEX_DIR:
    STAR_INDEX_DIR = os.path.join(
        ref_dir, f"STAR_{str(config['ref']['assembly']).strip()}"
    )
config["ref"]["star_index"] = STAR_INDEX_DIR


rule star_genome_generate:
    # Create STAR genome index for RNA-seq alignment
    input:
        fai   = config["ref"]["fasta"] + ".fai",
        gtf   = config["ref"]["gtf"],
        fasta = config["ref"]["fasta"]
    output:
        directory(STAR_INDEX_DIR)
    params:
        other_params = config["star_params"]["index"],
        mem_limit    = 85000000000
    conda:
        os.path.join(workflow.basedir, "envs", "star.yml")
    message:
        "Creating a STAR genome index"
    threads: 16
    resources:
        mem_mb = 90000  # ~85 GB for genome generation (matches limitGenomeGenerateRAM)
    log:
        os.path.join(STAR_INDEX_DIR, "genome_generate.log")
    benchmark:
        os.path.join(STAR_INDEX_DIR, "genome_generate.benchmark.txt")
    shell:
        """
        mkdir -p {output}
        mkdir -p $(dirname {log})

        # Calculate optimal genomeSAindexNbases based on genome size
        NUM_BASES=$(awk '{{sum += $2}} END {{
            val = (log(sum)/log(2))/2 - 1
            printf "%.0f", (val > 14 ? 14 : val)
        }}' {input.fai})

        STAR \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} \
            --genomeSAindexNbases $NUM_BASES \
            --limitGenomeGenerateRAM {params.mem_limit} \
            {params.other_params} > {log} 2>&1 || {{
                echo "[ERROR] STAR genome index generation failed." >> {log}
                exit 1
            }}
        """


rule star_align:
    # Perform RNA-seq alignment with STAR
    input:
        fq       = get_paired_trimmed_fq,
        star_idx = STAR_INDEX_DIR
    output:
        bam           = os.path.join("{outdir}", "bam", "{sample_id}.unsorted.bam"),
        tx_bam        = os.path.join("{outdir}", "bam", "{sample_id}.tx.bam"),
        log_final_out = os.path.join("{outdir}", "bam", "{sample_id}.Log.final.out")
    params:
        other_params  = lambda wc: config["star_params"].get(
            config["alignment"].get("param_type", "ffpe"), ""
        ),
        bam_sort_ram  = 32000000000,  # 32 GB for BAM sorting
        prefix        = lambda wc: os.path.join(wc.outdir, "bam", wc.sample_id + ".")
    conda:
        os.path.join(workflow.basedir, "envs", "star.yml")
    message:
        "{wildcards.sample_id}: Aligning with STAR"
    threads: 16
    resources:
        mem_mb = 40960  # STAR genome load + 32 GB BAM sort RAM
    log:
        os.path.join("{outdir}", "logs", "star", "{sample_id}.align.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.star_align.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {log})

        STAR \
            --genomeDir {input.star_idx} \
            --readFilesIn {input.fq[0]} {input.fq[1]} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.prefix} \
            --outSAMattrRGline ID:{wildcards.sample_id} SM:{wildcards.sample_id} \
            --readFilesCommand zcat \
            --limitBAMsortRAM {params.bam_sort_ram} \
            --outBAMcompression 6 \
            --runRNGseed 0 \
            {params.other_params} > {log} 2>&1 || {{
                echo "[ERROR] STAR alignment failed." >> {log}
                exit 1
            }}

        # Rename outputs; fail explicitly if any file is missing.
        # Log.final.out already matches the declared output path, so only validate it.
        mv {params.prefix}Aligned.out.bam                  {output.bam}    || {{ echo "[ERROR] Missing: Aligned.out.bam"                >> {log}; exit 1; }}
        mv {params.prefix}Aligned.toTranscriptome.out.bam  {output.tx_bam} || {{ echo "[ERROR] Missing: Aligned.toTranscriptome.out.bam" >> {log}; exit 1; }}
        [ -s {output.log_final_out} ] || {{ echo "[ERROR] Missing: Log.final.out" >> {log}; exit 1; }}
        mv {params.prefix}Log.out                          {params.prefix}Log.out.bak 2>/dev/null || true
        mv {params.prefix}Log.progress.out                 {params.prefix}Log.progress.out.bak 2>/dev/null || true
        mv {params.prefix}SJ.out.tab                       {params.prefix}SJ.out.tab.bak 2>/dev/null || true
        """


rule star_remove_shared_memory:
    # Remove STAR shared memory to free up system resources.
    # NOTE: This rule needs --genomeDir to target the correct genome;
    # without it, STAR removes whichever genome is loaded in shared memory.
    input:
        bam      = get_bam,
        star_idx = STAR_INDEX_DIR
    output:
        log_file = os.path.join("{outdir}", "bam", "{sample_id}.star_memory_removal.log")
    conda:
        os.path.join(workflow.basedir, "envs", "star.yml")
    message:
        "{wildcards.sample_id}: Removing STAR shared memory"
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.star_memory_removal.benchmark.txt")
    shell:
        """
        STAR --genomeDir {input.star_idx} --genomeLoad Remove > {output.log_file} 2>&1 || {{
            echo "[ERROR] Failed to remove STAR shared memory." >> {output.log_file}
            exit 1
        }}
        """


rule sort_bam:
    # Sort BAM file by coordinates and index.
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.unsorted.bam")
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.bam.bai")
    params:
        tempdir           = os.path.join("{outdir}", "bam"),
        memory_per_thread = "4G"
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "{wildcards.sample_id}: Sorting BAM by coordinates"
    threads: 10
    resources:
        mem_mb = 40960  # 4 GB per thread (matches memory_per_thread param)
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.sort.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.bam_sort.benchmark.txt")
    shell:
        """
        mkdir -p {params.tempdir}
        mkdir -p $(dirname {log})

        samtools sort \
            --write-index \
            -m {params.memory_per_thread} \
            -T {params.tempdir}/{wildcards.sample_id}.tmp \
            -@ {threads} \
            -o {output.bam}##idx##{output.bai} \
            {input.bam} > {log} 2>&1 || {{
                echo "[ERROR] BAM sorting failed." >> {log}
                exit 1
            }}

        if [ -s {output.bam} ] && [ -s {output.bai} ]; then
            echo "[INFO] Removing unsorted BAM to save space." >> {log}
            rm -f {input.bam}
        else
            echo "[ERROR] Sorted BAM or BAI missing." >> {log}
            exit 1
        fi
        """
