rule arriba:
    params:
        other_params = config['star_params']['fusion']
    input:
        get_paired_trimmed_fq,
        gtf = config['ref']['gtf'],
        fasta = config['ref']['fasta'],
        star_idx = config['ref']['star_idx']
    output:
        fusions = os.path.join("{outdir}", "arriba", "{sample_id}.fusions.tsv"),
        fusions_discarded = os.path.join("{outdir}", "arriba", "{sample_id}.fusions.discarded.tsv")
    message:
        "{wildcards.sample_id}: Calling fusion genes with Arriba."
    log:
        os.path.join("{outdir}", "arriba", "{sample_id}.arriba.log")
    conda:
        "../envs/arriba.yml"
    threads: 16
    shell:
        """
        STAR \
            --genomeDir {input.star_idx} \
            --outFileNamePrefix $(dirname {output.fusions})/{wildcards.sample_id}.fusionmapping. \
            --readFilesIn {input[0]} {input[1]} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --runRNGseed 0 \
            {params.other_params} 2>{log} \
            | arriba \
                -x /dev/stdin \
                -o {output.fusions} \
                -O {output.fusions_discarded} \
                -a {input.fasta} \
                -g {input.gtf} \
                -b $CONDA_PREFIX/var/lib/arriba/blacklist_{config[ref][assembly]}_*.tsv.gz \
                -k $CONDA_PREFIX/var/lib/arriba/known_fusions_{config[ref][assembly]}_*.tsv.gz \
                -t $CONDA_PREFIX/var/lib/arriba/known_fusions_{config[ref][assembly]}_*.tsv.gz \
                -p $CONDA_PREFIX/var/lib/arriba/protein_domains_{config[ref][assembly]}_*.gff3 >>{log} 2>&1
        """
