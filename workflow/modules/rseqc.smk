RSEQC_ENV = os.path.join(workflow.basedir, "envs", "rseqc.yml")
RSEQC_BED = get_ref_bed()


# Rule for BAM statistics
if is_rseqc_submodule_enabled("bam_stat"):
    rule bam_stat:
        input:
            bam = get_bam
        output:
            os.path.join("{outdir}", "rseqc", "bam_stat", "{sample_id}.bam_stat.txt")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.bam_stat.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 2048
        message:
            "{wildcards.sample_id}: Running BAM stat"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            bam_stat.py -i {input.bam} > {output} 2> {log} || {{ echo "[ERROR] bam_stat.py failed." >> {log}; exit 1; }}
            """

# Rule for inferring experiment (strand specificity)
if is_rseqc_submodule_enabled("infer_experiment"):
    rule infer_experiment:
        input:
            bam = get_bam,
            bed = RSEQC_BED
        output:
            os.path.join("{outdir}", "rseqc", "infer_experiment", "{sample_id}.infer_experiment.txt")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.infer_experiment.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 2048
        message:
            "{wildcards.sample_id}: Running Infer experiment"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log} || {{ echo "[ERROR] infer_experiment.py failed." >> {log}; exit 1; }}
            """

# Rule for inner distance analysis
if is_rseqc_submodule_enabled("inner_distance"):
    rule inner_distance:
        params:
           prefix = lambda w, output: output[0].replace(".inner_distance.txt", "")
        input:
            bam = get_bam,
            bed = RSEQC_BED
        output:
            os.path.join("{outdir}", "rseqc", "inner_distance", "{sample_id}.inner_distance.txt")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.inner_distance.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 2048
        message:
            "{wildcards.sample_id}: Running Inner distance"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            inner_distance.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1 || {{ echo "[ERROR] inner_distance.py failed." >> {log}; exit 1; }}
            """

# Rule for read distribution
if is_rseqc_submodule_enabled("read_distribution"):
    rule read_distribution:
        input:
            bam = get_bam,
            bed = RSEQC_BED
        output:
            os.path.join("{outdir}", "rseqc", "read_distribution", "{sample_id}.read_distribution.txt")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.read_distribution.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 2048
        message:
            "{wildcards.sample_id}: Running Read distribution"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log} || {{ echo "[ERROR] read_distribution.py failed." >> {log}; exit 1; }}
            """

# Rule for read duplication
if is_rseqc_submodule_enabled("read_duplication"):
    rule read_duplication:
        params:
            prefix = lambda wildcards: os.path.join(wildcards.outdir, "rseqc", "read_duplication", wildcards.sample_id)
        input:
            get_bam
        output:
            os.path.join("{outdir}", "rseqc", "read_duplication", "{sample_id}.DupRate_plot.pdf")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.read_duplication.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 2048
        message:
            "{wildcards.sample_id}: Running Read Duplication"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1 || {{ echo "[ERROR] read_duplication.py failed." >> {log}; exit 1; }}
            """

# Rule for read GC content
if is_rseqc_submodule_enabled("read_GC"):
    rule read_GC:
        params:
            prefix = lambda wildcards: os.path.join(wildcards.outdir, "rseqc", "read_GC", wildcards.sample_id)
        input:
            get_bam
        output:
            os.path.join("{outdir}", "rseqc", "read_GC", "{sample_id}.GC_plot.pdf")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.read_GC.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 2048
        message:
            "{wildcards.sample_id}: Running Read GC"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            read_GC.py -i {input} -o {params.prefix} > {log} 2>&1 || {{ echo "[ERROR] read_GC.py failed." >> {log}; exit 1; }}
            """

# Rule for junction annotation
if is_rseqc_submodule_enabled("junction_annotation"):
    rule junction_annotation:
        params:
            extra = "-q 255",  # STAR uses 255 as a score for unique mappers
            prefix = lambda wildcards: os.path.join(wildcards.outdir, "rseqc", "junction_annotation", wildcards.sample_id)
        input:
            bam = get_bam,
            bed = RSEQC_BED
        output:
            os.path.join("{outdir}", "rseqc", "junction_annotation", "{sample_id}.junction.bed")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.junction_annotation.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 4096
        message:
            "{wildcards.sample_id}: Running Junction annotation"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1 || {{ echo "[ERROR] junction_annotation.py failed." >> {log}; exit 1; }}
            """

# Rule for junction saturation
if is_rseqc_submodule_enabled("junction_saturation"):
    rule junction_saturation:
        params:
            extra = "-q 255",
            prefix = lambda wildcards: os.path.join(wildcards.outdir, "rseqc", "junction_saturation", wildcards.sample_id)
        input:
            bam = get_bam,
            bed = RSEQC_BED
        output:
            os.path.join("{outdir}", "rseqc", "junction_saturation", "{sample_id}.junctionSaturation_plot.pdf")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.junction_saturation.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 4096
        message:
            "{wildcards.sample_id}: Running Junction saturation"
        shell:
            """
            mkdir -p $(dirname {output})
            mkdir -p $(dirname {log})
            junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1 || {{ echo "[ERROR] junction_saturation.py failed." >> {log}; exit 1; }}
            """

# Rule for gene body coverage
if is_rseqc_submodule_enabled("gene_body_coverage"):
    rule gene_body_coverage:
        params:
            prefix = lambda wildcards: os.path.join(wildcards.outdir, "rseqc", "gene_body_coverage", wildcards.sample_id)
        input:
            bam = get_bam,
            bed = RSEQC_BED
        output:
            curves_pdf = os.path.join("{outdir}", "rseqc", "gene_body_coverage", "{sample_id}.geneBodyCoverage.curves.pdf"),
            txt = os.path.join("{outdir}", "rseqc", "gene_body_coverage", "{sample_id}.geneBodyCoverage.txt"),
            r_script = os.path.join("{outdir}", "rseqc", "gene_body_coverage", "{sample_id}.geneBodyCoverage.r")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.gene_body_coverage.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 4096
        message:
            "{wildcards.sample_id}: Running Gene body coverage"
        shell:
            """
            mkdir -p $(dirname {output.curves_pdf})
            mkdir -p $(dirname {log})
            geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.prefix} > /dev/null 2>&1 || {{ echo "[ERROR] geneBody_coverage.py failed." > {log}; exit 1; }}
            if [ -f "log.txt" ]; then
                mv log.txt {log}
            else
                echo "[WARNING] geneBody_coverage.py completed but log.txt was not generated." > {log}
            fi
            """

# Rule for transcript integrity number (TIN)
if is_rseqc_submodule_enabled("tin"):
    rule tin:
        params:
            outdir = lambda wildcards: os.path.join(wildcards.outdir, "rseqc", "tin")
        input:
            bam = get_bam,
            bai = get_bam_bai,
            bed = RSEQC_BED
        output:
            summary = os.path.join("{outdir}", "rseqc", "tin", "{sample_id}.tin.summary.txt"),
            tin_xls = os.path.join("{outdir}", "rseqc", "tin", "{sample_id}.tin.xls")
        log:
            os.path.join("{outdir}", "logs", "rseqc", "{sample_id}.tin.log")
        conda:
            RSEQC_ENV
        threads: 4
        resources:
            mem_mb = 4096
        message:
            "{wildcards.sample_id}: Running TIN (Transcript Integrity Number)"
        shell:
            """
            mkdir -p {params.outdir}
            mkdir -p $(dirname {log})
            BAM_ABS=$(readlink -f {input.bam})
            BED_ABS=$(readlink -f {input.bed})
            (
                cd {params.outdir}
                tin.py -i "$BAM_ABS" -r "$BED_ABS"
            ) > {log} 2>&1 || {{ echo "[ERROR] tin.py failed." >> {log}; exit 1; }}
            SUMMARY_FILE=$(ls {params.outdir}/*.summary.txt 2>/dev/null | head -n 1)
            TIN_FILE=$(ls {params.outdir}/*.tin.xls 2>/dev/null | head -n 1)
            [ -n "$SUMMARY_FILE" ] && mv -f "$SUMMARY_FILE" {output.summary} || (echo "[ERROR] TIN summary file not found." >> {log}; exit 1)
            [ -n "$TIN_FILE" ] && mv -f "$TIN_FILE" {output.tin_xls} || (echo "[ERROR] TIN xls file not found." >> {log}; exit 1)
            """
