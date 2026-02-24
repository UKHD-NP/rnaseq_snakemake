rule samtools_faidx:
    # Generate FASTA index
    input:
        config['ref']['fasta']
    output:
        config['ref']['fasta'] + ".fai"
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
    message:
        "Generating FASTA index"
    threads: 1
    shell:
        """
        samtools faidx {input}
        """

rule gffread:
    # Generate transcript FASTA from genome FASTA and GTF
    input:
        fasta = config['ref']['fasta'],
        gtf = config['ref']['gtf']
    output:
        config['ref']['tx_fasta']
    conda:
        os.path.join(workflow.basedir, "envs", "gffread.yml")
    message:
        "Generating transcript fasta"
    threads: 1
    log:
        os.path.join(ref_dir, "gffread", f"{config['ref']['assembly']}.log")
    shell:
        """
        echo "[INFO] Generating transcript FASTA with gffread..."
        echo "[INFO] Input FASTA: {input.fasta}"
        echo "[INFO] Input GTF: {input.gtf}"
        echo "[INFO] Output file: {output}"

        # Check if input files exist and have content
        if [ ! -s "{input.fasta}" ]; then
            echo "[ERROR] Input FASTA file is empty or does not exist: {input.fasta}"
            exit 1
        fi

        if [ ! -s "{input.gtf}" ]; then
            echo "[ERROR] Input GTF file is empty or does not exist: {input.gtf}"
            exit 1
        fi

        # Run with error checking
        gffread -g "{input.fasta}" -w "{output}" "{input.gtf}" 2> >(tee -a "{log}" >&2)

        # Check if output was created
        if [ ! -s "{output}" ]; then
            echo "[ERROR] gffread finished but output is missing or empty: {output}"
            exit 1
        fi
        """

rule gtf2bed:
    # Generate BED file from GTF
    input:
        config['ref']['gtf']
    output:
        config['ref']['bed']
    params:
        gtf2bed_script = os.path.join(workflow.basedir, "scripts", "gtf2bed")
    conda:
        os.path.join(workflow.basedir, "envs", "gtf2bed.yml")
    message:
        "Converting GTF to BED format"
    threads: 1
    log:
        os.path.join(ref_dir, "gtf2bed", f"{config['ref']['assembly']}.log")
    shell:
        """
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})

        if [ ! -f "{params.gtf2bed_script}" ]; then
            echo "[ERROR] gtf2bed script not found: {params.gtf2bed_script}" > {log}
            exit 1
        fi

        echo "[INFO] Using gtf2bed script: {params.gtf2bed_script}" > {log}
        perl "{params.gtf2bed_script}" "{input}" > "{output}" 2>> "{log}"

        # Check if output was created successfully
        if [ ! -s "{output}" ]; then
            echo "[ERROR] Failed to create BED file: {output}" >> {log}
            exit 1
        fi

        echo "[INFO] BED file created successfully: {output}" >> {log}
        """
