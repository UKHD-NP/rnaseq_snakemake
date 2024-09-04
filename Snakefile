import pandas as pd

# Load the config file
configfile: "config/config.yaml"

# Load the sample sheet from the config file
samplesheet = pd.read_csv(config["samples"])

# Access the trimming setting
trimming_active = config["trimming"]["enabled"]

# Create a new column that combines 'sample', 'unit', and 'condition'
samplesheet['sample_name'] = samplesheet['sample'] + "_" + samplesheet['unit'] + "_" + samplesheet['condition']

# Ensure that 'sample_name' column is properly recognized
sample_list = samplesheet['sample_name'].tolist()

def is_paired_end(sample_name):
    # Select the relevant row(s) for the given sample_name
    sample_units = samplesheet[samplesheet['sample_name'] == sample_name]
    
    # Check if 'fq2' is null (i.e., missing) for all the rows corresponding to this sample
    fq2_null = sample_units["fq2"].isnull()
    
    # Determine if this sample is paired-end
    paired = ~fq2_null
    
    # Check if all associated rows are either paired-end or single-end
    all_paired = paired.all()
    all_single = (~paired).all()
    
    # Assert that the sample is either completely paired-end or single-end
    assert (
        all_single or all_paired
    ), f"Invalid units for sample {sample_name}, must be all paired-end or all single-end."
    
    return all_paired

def get_paired_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.sample].iloc[0]
    if is_paired_end(wildcards.sample):
        return [sample_data['fq1'], sample_data['fq2']]  # Return a list of file paths
    else:
        raise ValueError(f"Sample {wildcards.sample} is not paired-end.")

def get_paired_trimmed_fq(wildcards):
    sample_data = samplesheet[samplesheet['sample_name'] == wildcards.sample].iloc[0]
    if is_paired_end(wildcards.sample):
        if config["trimming"]["enabled"]:
            return [
                f"results/trimmed/{wildcards.sample}/{wildcards.sample}_1.fastp.fastq.gz",
                f"results/trimmed/{wildcards.sample}/{wildcards.sample}_2.fastp.fastq.gz"
            ]
        else:
            return [sample_data['fq1'], sample_data['fq2']]
    else:
        raise ValueError(f"Sample {wildcards.sample} is not paired-end.")

# Add this rule to define the final target (all) for Snakemake
rule all:
    input:     
        "results/multiqc/multiqc_report.html"

# Rule for Preparing the Reference Genome
rule get_genome:
    input:
        config["ref"]["genome"]
    output:
        "resources/genome.fa"
    shell:
        "cp {input} {output}"

rule get_annotation:
    input:
        config["ref"]["annotation"]
    output:
        "resources/annotation.gtf"
    shell:
        "cp {input} {output}"

rule filter_gtf:
    input:
        gtf="resources/annotation.gtf",
        fasta="resources/genome.fa"
    output:
        gtf_filtered="resources/annotation.filtered.gtf"
    message: 
        "Refining the input GTF"
    log:
        "logs/filter_gtf.log"
    conda:
        "envs/python.yaml"
    shell:
        """
        python scripts/filter_gtf.py --gtf {input.gtf} --fasta {input.fasta} --prefix resources/annotation 2> {log}
        """

rule gtf_to_bed:
    input:
        gtf="resources/annotation.gtf"
    output:
        bed="resources/annotation.filtered.bed"
    message:
        "Converting a GTF file to a BED file format."
    log:
        "logs/gtf_to_bed.log"
    conda:
        "envs/perl.yaml"
    shell:
        "scripts/gtf2bed {input.gtf} > {output}"

rule rsem_prepare_reference:
    input:
        gtf="resources/annotation.filtered.gtf",
        fasta="resources/genome.fa"
    output:
        "resources/genome.transcripts.fa"
    message:
        "Preparing a reference transcript sequence for RNA-Seq quantification"
    params:
        prefix="rsem_genome"
    log:
        "logs/rsem_prepare_reference.log"
    conda:
        "envs/rsem.yaml"
    threads: 12
    shell:
        """
        rsem-prepare-reference \\
            --gtf {input.gtf} \\
            --num-threads {threads} \\
            {input.fasta} \\
            {params.prefix} 2> {log} &&
        cp {params.prefix}.transcripts.fa {output}
        """

rule samtools_faidx_sizes:
    input:
        fasta="resources/genome.fa"
    output:
        sizes="resources/genome.fa.sizes"
    message:
        "Listing the sizes of all sequences in the fasta file"
    log:
        "logs/samtools_faidx_sizes.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools faidx {input.fasta} 2> {log}
        cut -f 1,2 {input.fasta}.fai > {output}
        """

rule star_genome_generate:
    input:
        gtf="resources/annotation.filtered.gtf",
        fasta="resources/genome.fa"
    output:
        directory("resources/star_index")
    message:
        "Creating a reference genome index"
    log:
        "logs/star_genome_generate.log"
    conda:
        "envs/star.yaml"
    threads: 12
    shell:
        """
        samtools faidx {input.fasta} 2> {log}
        NUM_BASES=$(gawk '{{sum = sum + $2}}END{{if ((log(sum)/log(2))/2 - 1 > 14) {{printf "%.0f", 14}} else {{printf "%.0f", (log(sum)/log(2))/2 - 1}}}}' {input.fasta}.fai)
        STAR --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} --runThreadN {threads} --genomeSAindexNbases $NUM_BASES --limitGenomeGenerateRAM 77209411328 &>> {log}
        """

# Rule for Read QC
rule fastqc:
    input:
        get_paired_fq
    output:
        html1="results/fastqc/{sample}/{sample}_1.fastqc.html",
        zip1="results/fastqc/{sample}/{sample}_1.fastqc.zip",
        html2="results/fastqc/{sample}/{sample}_2.fastqc.html",
        zip2="results/fastqc/{sample}/{sample}_2.fastqc.zip"
    params:
        outdir="results/fastqc/{sample}",
        html1=lambda wildcards, output: output.html1.replace('__', '_'),
        zip1=lambda wildcards, output: output.zip1.replace('__', '_'),
        html2=lambda wildcards, output: output.html2.replace('__', '_'),
        zip2=lambda wildcards, output: output.zip2.replace('__', '_')
    message:
        "{wildcards.sample}: Performing quality control on raw reads"
    log:
        "logs/fastqc/{sample}.fastqc.log"
    conda:
        "envs/fastqc.yaml"    
    threads: 1
    resources:
        mem_mb = 1024
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run FastQC
        fastqc --quiet --threads {threads} {input[0]} {input[1]} \\
               --outdir {params.outdir} 2> {log}
        
        # Rename and move output files to match Snakemake expectations
        mv {params.outdir}/$(basename {input[0]} .fastq.gz)_fastqc.html  {output.html1}
        mv {params.outdir}/$(basename {input[0]} .fastq.gz)_fastqc.zip {output.zip1}
        mv {params.outdir}/$(basename {input[1]} .fastq.gz)_fastqc.html  {output.html2}
        mv {params.outdir}/$(basename {input[1]} .fastq.gz)_fastqc.zip {output.zip2}
        """

# Rule for Adapter and quality trimming
rule fastp:
    input:
        get_paired_fq
    output:
        out1="results/trimmed/{sample}/{sample}_1.fastp.fastq.gz",
        out2="results/trimmed/{sample}/{sample}_2.fastp.fastq.gz",
        json="results/trimmed/{sample}/{sample}.fastp.json",
        html="results/trimmed/{sample}/{sample}.fastp.html",
        unpaired1="results/trimmed/{sample}/{sample}_1.fail.fastq.gz",
        unpaired2="results/trimmed/{sample}/{sample}_2.fail.fastq.gz", 
    message:
        "{wildcards.sample}: Trimming and performing quality control on paired-end FASTQ files"
    log:
        "logs/fastp/{sample}.trimmed.log"
    conda:
        "envs/fastp.yaml"
    threads: 6
    shell: 
        """
        fastp \\
            --in1 {input[0]} \\
            --in2 {input[1]} \\
            --out1 {output.out1} \\
            --out2 {output.out2} \\
            --json {output.json} \\
            --html {output.html} \\
            --unpaired1 {output.unpaired1} \\
            --unpaired2 {output.unpaired2} \\
            --thread {threads} \\
            --detect_adapter_for_pe \\
            --trim_front1 1 --trim_front2 1 --length_required 50 \\
            2> {log}
        """

rule fastqc_after_trimming:
    input:
        fq1="results/trimmed/{sample}/{sample}_1.fastp.fastq.gz",
        fq2="results/trimmed/{sample}/{sample}_2.fastp.fastq.gz",
    output:
        html1="results/trimmed/{sample}/fastqc/{sample}_1.fastp_fastqc.html",
        zip1="results/trimmed/{sample}/fastqc/{sample}_1.fastp_fastqc.zip",
        html2="results/trimmed/{sample}/fastqc/{sample}_2.fastp_fastqc.html",
        zip2="results/trimmed/{sample}/fastqc/{sample}_2.fastp_fastqc.zip"
    params:
        outdir="results/trimmed/{sample}/fastqc",
    message:
        "{wildcards.sample}: Performing quality control on raw reads"
    log:
        "logs/fastqc/{sample}.fastqc.log"
    conda:
        "envs/fastqc.yaml"    
    threads: 1
    resources:
        mem_mb = 1024
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run FastQC
        fastqc --quiet --threads {threads} {input.fq1} {input.fq2} \\
               --outdir {params.outdir} 2> {log}
        """

# Rule for Alignment
rule align:
    input:
        unpack(get_paired_trimmed_fq),
        gtf="resources/annotation.filtered.gtf",
        star_index="resources/star_index"
    output:
        bam = "results/align/bam_original/{sample}/{sample}.Aligned.out.bam",
        transcriptome_bam = "results/align/bam_original/{sample}/{sample}.Aligned.toTranscriptome.out.bam",
        junction = "results/align/bam_original/{sample}/{sample}.SJ.out.tab"
    message:
        "{wildcards.sample}: Aligning paired-end FASTQ files"
    log:
        "logs/star_align/{sample}.star_align.log"
    conda:
        "envs/star.yaml"
    threads: 12
    params:
        outdir="results/align/bam_original/{sample}",
        star_common_params=config["star"]["common_params"]
    shell:
        """
        # Ensure output directory exists
        mkdir -p {params.outdir}

        # Run STAR alignment
        STAR --genomeDir {input.star_index} \\
             --readFilesIn {input[0]} {input[1]} \\
             --runThreadN {threads} \\
             --genomeLoad NoSharedMemory \\
             --outFileNamePrefix {params.outdir}/{wildcards.sample}. \\
             --sjdbGTFfile {input.gtf} \\
             --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \\
             --twopassMode Basic \\
             --outSAMtype BAM Unsorted \\
             --readFilesCommand zcat \\
             --runRNGseed 0 \\
             --quantMode TranscriptomeSAM \\
             --outSAMstrandField intronMotif \\
             --outSAMattributes NH HI NM MD AS XS \\
             {params.star_common_params} 2> {log}
        """
# Rule for Sort and index alignments 
rule samtools_process:
    input:
        bam="results/align/bam_original/{sample}/{sample}.Aligned.out.bam",
        fasta="resources/genome.fa"
    output:
        sorted_bam="results/align/bam_original/{sample}/{sample}.sorted.bam",
        sorted_bam_bai="results/align/bam_original/{sample}/{sample}.sorted.bam.bai",
        stats="results/align/bam_original/{sample}/{sample}.sorted.bam.stats",
        flagstat="results/align/bam_original/{sample}/{sample}.sorted.bam.flagstat",
        idxstats="results/align/bam_original/{sample}/{sample}.sorted.bam.idxstats"
    message:
        "{wildcards.sample}: Sorting and indexing alignments"
    log:
        "logs/samtools_process/{sample}.samtools_process.log"
    conda:
        "envs/samtools.yaml"
    threads: 8
    shell:
        """
        # Sorts alignments in BAM file
        samtools sort -@ {threads} -m 8G -T tmp -o {output.sorted_bam} {input.bam} 2> {log}
        
        # Indexes the sorted BAM file
        samtools index -@ 1 {output.sorted_bam} 2>> {log}
        
        # Generates statistics for the sorted BAM file
        samtools stats --threads 1 --reference {input.fasta} {output.sorted_bam} > {output.stats} 2>> {log}
        
        # Counts the number of alignments in the BAM file for each FLAG type
        samtools flagstat --threads 1 {output.sorted_bam} > {output.flagstat} 2>> {log}
        
        # Produces index statistics for the sorted BAM file
        samtools idxstats {output.sorted_bam} > {output.idxstats} 2>> {log}
        """

# Rule for Duplicate read marking
rule picard_mark_duplicates:
    input:
        bam="results/align/bam_original/{sample}/{sample}.sorted.bam",
        fasta="resources/genome.fa"
    output:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        metrics="results/align/bam_markup/{sample}/{sample}.markdup.sorted.MarkDuplicates.metrics.txt"
    message:
        "{wildcards.sample}: Marking duplicate reads"
    log:
        "logs/markdup/{sample}.markduplicates.log"
    conda:
        "envs/picard_markduplicates.yaml"
    resources:
        mem_mb=4096
    threads: 12
    shell:
        """
        picard MarkDuplicates \\
        --SORTING_COLLECTION_SIZE_RATIO 0.12 \\
        --ASSUME_SORTED true \\
        --REMOVE_DUPLICATES false \\
        --VALIDATION_STRINGENCY LENIENT \\
        --TMP_DIR tmp \\
        --INPUT {input.bam} \\
        --OUTPUT {output.bam} \\
        --REFERENCE_SEQUENCE {input.fasta} \\
        --METRICS_FILE {output.metrics} \\
        2> {log}
        """

rule samtools_markdup_processing:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        fasta="resources/genome.fa"
    output:
        bam_bai="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.bai",
        bam_stats="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.stats",
        bam_flagstat="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.flagstat",
        bam_idxstats="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.idxstats"
    message:
        "{wildcards.sample}: Running Samtools to mark duplicate reads"
    log:
        "logs/markdup/{sample}.markduplicates.samtools.log"
    conda:
        "envs/samtools.yaml"
    threads: 1
    shell:
        """
        # Index the BAM file
        samtools index -@ {threads} {input.bam} {output.bam_bai} 2> {log}
        
        # Get BAM statistics
        samtools stats --threads {threads} -r {input.fasta} {input.bam} > {output.bam_stats} 2>> {log}
        
        # Get BAM flagstat
        samtools flagstat --threads {threads} {input.bam} > {output.bam_flagstat} 2>> {log}
        
        # Get BAM idxstats
        samtools idxstats {input.bam} > {output.bam_idxstats} 2>> {log}
        """

# Rule for Transcript assembly and quantification
rule stringtie:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf"
    output:
        gtf="results/stringtie/{sample}/{sample}.transcripts.gtf",
        abundance="results/stringtie/{sample}/{sample}.gene.abundance.txt",
        coverage="results/stringtie/{sample}/{sample}.coverage.gtf",
        ballgown_dir=directory("results/stringtie/{sample}/{sample}.ballgown")
    message:
        "{wildcards.sample}: Running Stringtie to assemble and quantify transcripts"
    log:
        "logs/stringtie/{sample}.stringtie.log"
    conda:
        "envs/stringtie.yaml"
    threads: 6
    params:
        stringtie_enabled=config.get("stringtie", {}).get("enabled", True)
    shell:
        """
        if [ "{params.stringtie_enabled}" = "True" ]; then
            stringtie {input.bam} \\
                      --rf \\
                      -G {input.gtf} \\
                      -o {output.gtf} \\
                      -A {output.abundance} \\
                      -C {output.coverage} \\
                      -b {output.ballgown_dir} \\
                      -p {threads} \\
                      -v \\
                      -e 2> {log}
        else
            mkdir -p $(dirname {output.gtf})
            mkdir -p $(dirname {output.abundance})
            mkdir -p $(dirname {output.coverage})
            mkdir -p {output.ballgown_dir}
            touch {output.gtf}
            touch {output.abundance}
            touch {output.coverage}
            echo "Stringtie disabled in configuration. Skipping..." > {log}
        fi
        """

# Rule for Read counting relative to gene biotype
rule feature_counts:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf",
        header="templates/biotypes_header.txt"
    output:
        counts="results/feature_counts/{sample}/{sample}.featurecounts.txt",
        summary="results/feature_counts/{sample}/{sample}.featurecounts.txt.summary",
        stats="results/feature_counts/{sample}/{sample}.biotype_counts_rrna_mqc.tsv"
    message:
        "{wildcards.sample}: Read counting relative to gene biotype"
    log:
        count_log="logs/feature_counts/{sample}.feature_counts.log",
        stats_log="logs/feature_counts/{sample}.feature_counts.stats.log"
    conda:
        "envs/featurecounts.yaml"
    threads: 6
    params:
        feature_counts_enabled=config.get("feature_counts", {}).get("enabled", True),
        featurecounts_group_type=config["feature_counts"].get("group_type", "gene_id")
    shell:
        """
        if [ "{params.feature_counts_enabled}" = "True" ]; then
            featureCounts -B -C \\
                -g {params.featurecounts_group_type} \\
                -t exon \\
                -p \\
                -T {threads} \\
                -a {input.gtf} \\
                -s 2 \\
                -o {output.counts} \\
                {input.bam} &> {log.count_log} &&

            cut -f 1,7 {output.counts} \\
                | tail -n +3 \\
                | cat {input.header} - \\
                > {output.stats} &&

            python scripts/mqc_features_stat.py \\
                {output.stats} \\
                -s {wildcards.sample} \\
                -f rRNA \\
                -o {output.stats} &> {log.stats_log}
        else
            mkdir -p $(dirname {output.counts})
            mkdir -p $(dirname {output.summary})
            mkdir -p $(dirname {output.stats})
            touch {output.counts}
            touch {output.summary}
            touch {output.stats}
            echo "FeatureCounts disabled in configuration. Skipping..." > {log.count_log}
        fi
        """


# Rule for Qualimap
rule qualimap:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf",
    output:
        "results/qualimap/{sample}/qualimapReport.html"    
    message:
        "{wildcards.sample}: Running Qualimap"
    log:
        "logs/qualimap/{sample}.qualimap.log"
    conda:
        "envs/qualimap.yaml"
    threads: 6
    params:
        qualimap_enabled=config.get("qualimap", {}).get("enabled", True)
    shell:
        """
        if [ "{params.qualimap_enabled}" = "True" ]; then
            unset DISPLAY &&
            mkdir -p tmp &&
            export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp &&
            qualimap --java-mem-size=4G rnaseq \\
            -bam {input.bam} \\
            -gtf {input.gtf} \\
            -p strand-specific-reverse \\
            -pe  \\
            -outdir results/qualimap/{wildcards.sample} 2> {log}
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Qualimap disabled in configuration. Skipping..." > {log}
        fi
        """

# Rule for Assessment of technical / biological read duplication
rule dupradar:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        gtf="resources/annotation.filtered.gtf"
    output:
        "results/dupradar/{sample}/{sample}_duprateExpDens.pdf"
    log:
        "logs/dupradar/{sample}.dupradar.log"
    params:
        prefix=lambda wildcards, output: output[0].replace("_duprateExpDens.pdf", ""),
        dupradar_enabled=config.get("dupradar", {}).get("enabled", True)
    message:
        "{wildcards.sample}: Running dupRadar to evaluate technical and biological read duplication"
    conda:
        "envs/dupradar.yaml"
    threads: 4
    shell:
        """
        if [ "{params.dupradar_enabled}" = "True" ]; then
            scripts/dupradar.r {input.bam} {params.prefix} {input.gtf} 2 paired {threads} 2> {log}
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "dupRadar disabled in configuration. Skipping..." > {log}
        fi
        """

# Rule for RSeQC 
rule bam_stat:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam"
    output:
        "results/rseqc/bam_stat/{sample}/{sample}.bam_stat.txt"
    log:
        "logs/rseqc/{sample}.bam_stat.log"
    conda:
        "envs/rseqc.yaml"
    threads: 1
    params:
        enabled=config.get("rseqc", {}).get("bam_stat", {}).get("enabled", True)
    message: 
        "{wildcards.sample}: Running BAM stat"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            bam_stat.py -i {input.bam} > {output} 2> {log}
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "BAM stat disabled in configuration. Skipping..." > {log}
        fi
        """

rule infer_experiment:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        "results/rseqc/infer_experiment/{sample}/{sample}.infer_experiment.txt"
    log:
        "logs/rseqc/{sample}.infer_experiment.log"
    conda:
        "envs/rseqc.yaml"
    threads: 1
    params:
        enabled=config.get("rseqc", {}).get("infer_experiment", {}).get("enabled", True)
    message: 
        "{wildcards.sample}: Running Infer experiment"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Infer experiment disabled in configuration. Skipping..." > {log}
        fi
        """

rule inner_distance:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        "results/rseqc/inner_distance/{sample}/{sample}.inner_distance.txt"
    priority: 1
    log:
        "logs/rseqc/{sample}.inner_distance.log"
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
        enabled=config.get("rseqc", {}).get("inner_distance", {}).get("enabled", True)
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running Inner distance"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            inner_distance.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Inner distance disabled in configuration. Skipping..." > {log}
        fi
        """

rule read_distribution:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        "results/rseqc/read_distribution/{sample}/{sample}.read_distribution.txt"
    priority: 1
    log:
        "logs/rseqc/{sample}.read_distribution.log"
    params:
        enabled=config.get("rseqc", {}).get("read_distribution", {}).get("enabled", True)
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running Read distribution"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Read distribution disabled in configuration. Skipping..." > {log}
        fi
        """

rule read_duplication:
    input:
        "results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam"
    output:
        "results/rseqc/read_duplication/{sample}/{sample}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        "logs/rseqc/{sample}.read_duplication.log"
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
        enabled=config.get("rseqc", {}).get("read_duplication", {}).get("enabled", True)
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running Read Duplication"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Read Duplication disabled in configuration. Skipping..." > {log}
        fi
        """

rule read_GC:
    input:
        "results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam"
    output:
        "results/rseqc/read_GC/{sample}/{sample}.readgc.GC_plot.pdf"
    priority: 1
    log:
        "logs/rseqc/{sample}.read_GC.log"
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
        enabled=config.get("rseqc", {}).get("read_GC", {}).get("enabled", True)
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running Read GC"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            read_GC.py -i {input} -o {params.prefix} > {log} 2>&1
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Read GC disabled in configuration. Skipping..." > {log}
        fi
        """

rule junnction_annotation:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed"
    output:
        "results/rseqc/junction_annotation/{sample}/{sample}.junction.bed"
    priority: 1
    log:
        "logs/rseqc/{sample}.junction_annotation.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
        enabled=config.get("rseqc", {}).get("junnction_annotation", {}).get("enabled", True)
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running Junction annotation"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Junction annotation disabled in configuration. Skipping..." > {log}
        fi
        """

rule junction_saturation:
    input:
        bam="results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam",
        bed="resources/annotation.filtered.bed",
    output:
        "results/rseqc/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf"
    priority: 1
    log:
        "logs/rseqc/{sample}.junction_saturation.log"
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
        enabled=config.get("rseqc", {}).get("junction_saturation", {}).get("enabled", True)
    conda:
        "envs/rseqc.yaml"
    threads: 1
    message: 
        "{wildcards.sample}: Running Junction saturation"
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1
        else
            mkdir -p $(dirname {output})
            touch {output}
            echo "Junction saturation disabled in configuration. Skipping..." > {log}
        fi
        """

rule multiqc:
    input:
        expand("results/fastqc/{sample}/{sample}_1.fastqc.html", sample=sample_list),
        expand("results/fastqc/{sample}/{sample}_1.fastqc.zip", sample=sample_list),
        expand("results/fastqc/{sample}/{sample}_2.fastqc.html", sample=sample_list),
        expand("results/fastqc/{sample}/{sample}_2.fastqc.zip", sample=sample_list),
        expand("results/trimmed/{sample}/{sample}.fastp.json", sample=sample_list),
        expand("results/trimmed/{sample}/{sample}.fastp.html",sample=sample_list),
        expand("results/trimmed/{sample}/fastqc/{sample}_1.fastp_fastqc.html", sample=sample_list),
        expand("results/trimmed/{sample}/fastqc/{sample}_1.fastp_fastqc.zip", sample=sample_list),
        expand("results/trimmed/{sample}/fastqc/{sample}_2.fastp_fastqc.html", sample=sample_list),
        expand("results/trimmed/{sample}/fastqc/{sample}_2.fastp_fastqc.zip", sample=sample_list),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam.stats", sample=sample_list),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam.flagstat", sample=sample_list),
        expand("results/align/bam_original/{sample}/{sample}.sorted.bam.idxstats", sample=sample_list),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.stats", sample=sample_list),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.flagstat", sample=sample_list),
        expand("results/align/bam_markup/{sample}/{sample}.markdup.sorted.bam.idxstats", sample=sample_list),
        expand("results/stringtie/{sample}/{sample}.gene.abundance.txt", sample=sample_list),
        expand("results/feature_counts/{sample}/{sample}.featurecounts.txt.summary", sample=sample_list),
        expand("results/feature_counts/{sample}/{sample}.biotype_counts_rrna_mqc.tsv", sample=sample_list),
        expand("results/dupradar/{sample}/{sample}_duprateExpDens.pdf", sample=sample_list),
        expand("results/qualimap/{sample}/qualimapReport.html", sample=sample_list),
        expand("results/rseqc/bam_stat/{sample}/{sample}.bam_stat.txt", sample=sample_list),
        expand("results/rseqc/infer_experiment/{sample}/{sample}.infer_experiment.txt", sample=sample_list),
        expand("results/rseqc/inner_distance/{sample}/{sample}.inner_distance.txt", sample=sample_list),
        expand("results/rseqc/read_distribution/{sample}/{sample}.read_distribution.txt", sample=sample_list),
        expand("results/rseqc/read_duplication/{sample}/{sample}.readdup.DupRate_plot.pdf", sample=sample_list),
        expand("results/rseqc/read_GC/{sample}/{sample}.readgc.GC_plot.pdf", sample=sample_list),
        expand("results/rseqc/junction_saturation/{sample}/{sample}.junctionSaturation_plot.pdf", sample=sample_list),
        multiqc_config="templates/multiqc_config.yml"
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    threads: 1
    message: 
        "Running MultiQC on all samples."
    shell:
        "multiqc {input} --outdir results/multiqc --config {input.multiqc_config}"





