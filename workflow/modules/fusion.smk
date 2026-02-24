# Arriba fusion detection
if is_enabled("fusion"):
    rule arriba:
        params:
            other_params = config['star_params']['fusion']
        input:
            fastq = get_paired_trimmed_fq,
            gtf = config['ref']['gtf'],
            fasta = config['ref']['fasta'],
            star_idx = config['ref']['star_index']
        output:
            fusions = os.path.join("{outdir}", "arriba", "{sample_id}.fusions.tsv"),
            fusions_discarded = os.path.join("{outdir}", "arriba", "{sample_id}.fusions.discarded.tsv")
        message:
            "{wildcards.sample_id}: Calling fusion genes with Arriba."
        log:
            os.path.join("{outdir}", "arriba", "{sample_id}.arriba.log")
        conda:
            os.path.join(workflow.basedir, "envs", "arriba.yml")
        threads: 16
        resources:
            mem_mb = 50000  # STAR + Arriba: memory-intensive fusion detection
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
                    -b $CONDA_PREFIX/var/lib/arriba/blacklist_{config[ref][assembly]}_*.tsv.gz \
                    -k $CONDA_PREFIX/var/lib/arriba/known_fusions_{config[ref][assembly]}_*.tsv.gz \
                    -t $CONDA_PREFIX/var/lib/arriba/known_fusions_{config[ref][assembly]}_*.tsv.gz \
                    -p $CONDA_PREFIX/var/lib/arriba/protein_domains_{config[ref][assembly]}_*.gff3 >>{log} 2>&1 || {{ echo "[ERROR] Arriba fusion detection failed." >> {log}; exit 1; }}

            if [ ! -s "{output.fusions}" ] || [ ! -s "{output.fusions_discarded}" ]; then
                echo "[ERROR] Arriba outputs are missing or empty." >> {log}
                exit 1
            fi
            """
