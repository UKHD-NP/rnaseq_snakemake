# Arriba fusion detection
rule arriba:
    input:
        fastq = get_paired_align_fq,
        gtf = config['ref']['gtf'],
        fasta = config['ref']['fasta'],
        star_idx = config['ref']['star_index']
    output:
        fusions = os.path.join("{outdir}", "arriba", "{sample_id}.fusions.tsv"),
        fusions_discarded = os.path.join("{outdir}", "arriba", "{sample_id}.fusions.discarded.tsv")
    params:
        other_params = config['star_params']['fusion']
    conda:
        os.path.join(workflow.basedir, "envs", "arriba.yml")
    message:
        "{wildcards.sample_id}: Calling fusion genes with Arriba."
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: min(130000 + (attempt - 1) * 35000, 200000),
        runtime = lambda wildcards, attempt: attempt * 480
    log:
        os.path.join("{outdir}", "arriba", "{sample_id}.arriba.log")
    shell:
        """
        set -o pipefail
        mkdir -p $(dirname {output.fusions})
        mkdir -p $(dirname {log})

        STAR \
            --genomeDir {input.star_idx} \
            --outFileNamePrefix $(dirname {output.fusions})/{wildcards.sample_id}.fusionmapping. \
            --readFilesIn {input.fastq[0]} {input.fastq[1]} \
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
                -b $CONDA_PREFIX/var/lib/arriba/blacklist_{config[fusion][assembly]}_*.tsv.gz \
                -k $CONDA_PREFIX/var/lib/arriba/known_fusions_{config[fusion][assembly]}_*.tsv.gz \
                -t $CONDA_PREFIX/var/lib/arriba/known_fusions_{config[fusion][assembly]}_*.tsv.gz \
                -p $CONDA_PREFIX/var/lib/arriba/protein_domains_{config[fusion][assembly]}_*.gff3 >>{log} 2>&1 || {{ echo "[ERROR] Arriba fusion detection failed." >> {log}; exit 1; }}

        if [ ! -s "{output.fusions}" ] || [ ! -s "{output.fusions_discarded}" ]; then
            echo "[ERROR] Arriba outputs are missing or empty." >> {log}
            exit 1
        fi
        """
