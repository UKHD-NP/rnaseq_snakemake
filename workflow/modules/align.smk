rule star_genome_generate:
    input:
        gtf = config['ref']['gtf'],
        fasta = config['ref']['fasta']
    output:
        directory(config['ref']['star_idx'])
    message:
        "Creating a STAR genome index"
    log:
        os.path.join(config['ref']['star_idx'], "genome_generate.log")
    conda:
        "../envs/star.yaml"
    threads: 12
    shell:
        """
        samtools faidx {input.fasta} 2> {log}
        NUM_BASES=$(gawk '{{sum = sum + $2}}END{{if ((log(sum)/log(2))/2 - 1 > 14) {{printf "%.0f", 14}} else {{printf "%.0f", (log(sum)/log(2))/2 - 1}}}}' {input.fasta}.fai)
        STAR \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} \
            --genomeSAindexNbases $NUM_BASES &>> {log}
        """


rule align:
    params:
        other_params = config['star_params']['ffpe']
    input:
        get_paired_trimmed_fq,
        star_idx = config['ref']['star_idx']
    output:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.bam.bai"),
    message:
        "{wildcards.sample_id}: Aligning with STAR"
    conda:
        "../envs/star.yaml"
    threads: 12
    shell:
        """
        STAR --genomeDir {input.star_idx} \
             --readFilesIn {input[0]} {input[1]} \
             --runThreadN {threads} \
             --outFileNamePrefix {wildcards.outdir}/bam/ \
             --outSAMattrRGline ID:{wildcards.sample_id} SM:{wildcards.sample_id} \
             --readFilesCommand zcat \
             --runRNGseed 0 \
             {params.other_params} > /dev/null

             mv {wildcards.outdir}/bam/Aligned.sortedByCoord.out.bam {output.bam}
             samtools index {output.bam}
        
        mv {wildcards.outdir}/bam/Log.out {wildcards.outdir}/bam/{wildcards.sample_id}.Log.out
        mv {wildcards.outdir}/bam/Log.progress.out {wildcards.outdir}/bam/{wildcards.sample_id}.Log.progress.out
        mv {wildcards.outdir}/bam/Log.final.out {wildcards.outdir}/bam/{wildcards.sample_id}.Log.final.out
        mv {wildcards.outdir}/bam/SJ.out.tab {wildcards.outdir}/bam/{wildcards.sample_id}.SJ.out.tab
        """
