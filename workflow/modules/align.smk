rule star_genome_generate:
    # Create STAR genome index for RNA-seq alignment
    params:
        other_params = config['star_params']['index'],
        mem_limit = 85000000000
    input:
        fai = config['ref']['fasta'] + ".fai",
        gtf = config['ref']['gtf'],
        fasta = config['ref']['fasta']
    output:
        directory(config['ref']['star_idx'])
    message:
        "Creating a STAR genome index"
    log:
        os.path.join(config['ref']['star_idx'], "genome_generate.log")
    conda:
        os.path.join(workflow.basedir, "envs", "star.yml")
    threads: 16
    benchmark:
        os.path.join(config['ref']['star_idx'], "genome_generate.benchmark.txt")
    shell:
        """
        # Calculate optimal genomeSAindexNbases based on genome size
        mkdir -p {output}
        mkdir -p $(dirname {log})
        NUM_BASES=$(gawk '{{sum = sum + $2}}END{{if ((log(sum)/log(2))/2 - 1 > 14) {{printf "%.0f", 14}} else {{printf "%.0f", (log(sum)/log(2))/2 - 1}}}}' {input.fai})
        
        STAR \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} \
            --genomeSAindexNbases $NUM_BASES \
            --limitGenomeGenerateRAM {params.mem_limit} \
            {params.other_params} &> {log} || {{ echo "[ERROR] STAR genome index generation failed." >> {log}; exit 1; }}
        """

rule star_align:
    # Perform RNA-seq alignment with STAR
    params:
        other_params = lambda wildcards: config['star_params'].get(
            config['alignment'].get('param_type', 'ffpe'),
            ''
        ),
        bam_sort_ram = 32000000000  # 32GB for BAM sorting
    input:
        get_paired_trimmed_fq,
        star_idx = config['ref']['star_idx']
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.unsorted.bam"),
        tx_bam = os.path.join("{outdir}", "bam", "{sample_id}.tx.bam"),
        log_final_out = os.path.join("{outdir}", "bam", "{sample_id}.Log.final.out")
    message:
        "{wildcards.sample_id}: Aligning with STAR"
    conda:
        os.path.join(workflow.basedir, "envs", "star.yml")
    threads: 16
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.star_align.benchmark.txt")
    log:
        os.path.join("{outdir}", "logs", "star", "{sample_id}.align.log")
    shell:
        """
        mkdir -p {wildcards.outdir}/bam
        mkdir -p $(dirname {log})

        STAR --genomeDir {input.star_idx} \
             --readFilesIn {input[0]} {input[1]} \
             --runThreadN {threads} \
             --outFileNamePrefix {wildcards.outdir}/bam/ \
             --outSAMattrRGline ID:{wildcards.sample_id} SM:{wildcards.sample_id} \
             --readFilesCommand zcat \
             --limitBAMsortRAM {params.bam_sort_ram} \
             --outBAMcompression 6 \
             --runRNGseed 0 \
             {params.other_params} &> {log} || {{ echo "[ERROR] STAR alignment failed." >> {log}; exit 1; }}

        # Rename output files for consistency
        mv {wildcards.outdir}/bam/Aligned.out.bam {output.bam} || {{ echo "[ERROR] Missing STAR output: Aligned.out.bam." >> {log}; exit 1; }}
        mv {wildcards.outdir}/bam/Log.out {wildcards.outdir}/bam/{wildcards.sample_id}.Log.out || {{ echo "[ERROR] Missing STAR output: Log.out." >> {log}; exit 1; }}
        mv {wildcards.outdir}/bam/Log.progress.out {wildcards.outdir}/bam/{wildcards.sample_id}.Log.progress.out || {{ echo "[ERROR] Missing STAR output: Log.progress.out." >> {log}; exit 1; }}
        mv {wildcards.outdir}/bam/Log.final.out {output.log_final_out} || {{ echo "[ERROR] Missing STAR output: Log.final.out." >> {log}; exit 1; }}
        mv {wildcards.outdir}/bam/SJ.out.tab {wildcards.outdir}/bam/{wildcards.sample_id}.SJ.out.tab || {{ echo "[ERROR] Missing STAR output: SJ.out.tab." >> {log}; exit 1; }}
        mv {wildcards.outdir}/bam/Aligned.toTranscriptome.out.bam {output.tx_bam} || {{ echo "[ERROR] Missing STAR output: Aligned.toTranscriptome.out.bam." >> {log}; exit 1; }}
        """

rule star_remove_shared_memory:
    # Remove STAR shared memory to free up system resources
    input:
        bam = get_bam
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
        # Run STAR and directly output to the log file
        STAR --genomeLoad Remove &> {output.log_file} || {{ 
            echo "[ERROR] Failed to remove STAR shared memory." >> {output.log_file}
            exit 1
        }}
        """

rule sort_bam:
    # Sort BAM file by coordinates and index
    params:
        tempdir = os.path.join("{outdir}", "bam"),
        memory_per_thread = "4G"
    input:
        os.path.join("{outdir}", "bam", "{sample_id}.unsorted.bam")
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.bam.bai")
    message:
        "{wildcards.sample_id}: Sorting BAM by Coordinates"
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.sort.log")
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    threads: 10
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
            {input} >{log} 2>&1

        if [ -s "{output.bam}" ]; then
            echo "[INFO] Removing input BAM to save space." >> {log}
            rm -f {input}
        else
            echo "[ERROR] BAM sorting failed." >> {log}
            exit 1
        fi
        """
        
rule samtools_stats:
    # Generate comprehensive statistics for BAM files
    # Works with both markdup BAMs (if enabled) and regular sorted BAMs (if markdup disabled)
    input:
        bam = get_bam,
        fasta = config['ref']['fasta']
    output:
        bam_stats = os.path.join("{outdir}", "bam", "{sample_id}.bam.stats"),
        bam_flagstat = os.path.join("{outdir}", "bam", "{sample_id}.bam.flagstat"),
        bam_idxstats = os.path.join("{outdir}", "bam", "{sample_id}.bam.idxstats")
    message:
        "{wildcards.sample_id}: Running Samtools statistics"
    log:
        os.path.join("{outdir}", "logs", "samtools", "{sample_id}.samtools_stats.log")
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    threads: 4
    benchmark:
        os.path.join("{outdir}", "benchmarks", "{sample_id}.samtools_stats.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.bam_stats})
        mkdir -p $(dirname {log})

        # Generate comprehensive BAM statistics
        samtools stats --threads {threads} -r {input.fasta} {input.bam} > {output.bam_stats} 2>> {log} || {{ echo "[ERROR] samtools stats failed." >> {log}; exit 1; }}

        # Create flagstat summary
        samtools flagstat --threads 2 {input.bam} > {output.bam_flagstat} 2>> {log} || {{ echo "[ERROR] samtools flagstat failed." >> {log}; exit 1; }}

        # Generate chromosome-level read mapping statistics
        samtools idxstats {input.bam} > {output.bam_idxstats} 2>> {log} || {{ echo "[ERROR] samtools idxstats failed." >> {log}; exit 1; }}
        """
