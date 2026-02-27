# Rule for Transcript assembly and quantification
rule stringtie:
    input:
        bam = get_bam,
        gtf = config['ref']['gtf']
    output:
        gtf = os.path.join("{outdir}", "stringtie", "{sample_id}.transcripts.gtf"),
        abundance = os.path.join("{outdir}", "stringtie", "{sample_id}.gene.abundance.txt"),
        coverage = os.path.join("{outdir}", "stringtie", "{sample_id}.coverage.gtf"),
        ballgown_dir = directory(os.path.join("{outdir}", "stringtie", "{sample_id}.ballgown"))
    params:
        strand_flag = {
            "rf": "--rf",
            "fr": "--fr",
        }.get(str(config.get("stringtie", {}).get("strand", "rf")).lower(), "")
    conda:
        os.path.join(workflow.basedir, "envs", "stringtie.yml")
    message:
        "{wildcards.sample_id}: Running Stringtie to assemble and quantify transcripts"
    threads: 6
    resources:
        mem_mb = 8192
    log:
        os.path.join("{outdir}", "logs", "stringtie", "{sample_id}.stringtie.log")
    benchmark:
        os.path.join("{outdir}", "benchmarks", "stringtie.{sample_id}.benchmark.txt")
    shell:
        """
        mkdir -p $(dirname {output.gtf})
        mkdir -p {output.ballgown_dir}
        mkdir -p $(dirname {log})

        stringtie {input.bam} \
            {params.strand_flag} \
            -G {input.gtf} \
            -o {output.gtf} \
            -A {output.abundance} \
            -C {output.coverage} \
            -b {output.ballgown_dir} \
            -p {threads} \
            -v \
            -e 2> {log} || {{ echo "[ERROR] StringTie failed." >> {log}; exit 1; }}
        """
