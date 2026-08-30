rule delete_tmp:
    # Clean up temporary files
    input:
        bam = os.path.join("{outdir}", "bam", "{sample_id}.bam"),
        bai = os.path.join("{outdir}", "bam", "{sample_id}.bam.bai"),
        fastqc = lambda wildcards: [
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_1_fastqc.zip"),
            os.path.join(wildcards.outdir, "trim", f"{wildcards.sample_id}_trimmed_2_fastqc.zip"),
        ] if is_enabled("trimming") else []
    output:
        log = os.path.join("{outdir}", "logs", "{sample_id}.deletion.log")
    params:
        bam_dir = os.path.join("{outdir}", "bam"),
        fq1 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_1.fastq.gz"),
        fq2 = os.path.join("{outdir}", "trim", "{sample_id}_trimmed_2.fastq.gz"),
        fq1_fail = os.path.join("{outdir}", "trim", "{sample_id}_1.fail.fastq.gz"),
        fq2_fail = os.path.join("{outdir}", "trim", "{sample_id}_2.fail.fastq.gz"),
        raw_fq1 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_1.fastq.gz"),
        raw_fq2 = os.path.join("{outdir}", "raw_merged", "{sample_id}_merged_2.fastq.gz"),
        raw_dir = os.path.join("{outdir}", "raw_merged"),
        delete_trimming = lambda wildcards: str(
            is_enabled("trimming") and
            as_bool(config.get("trimming", {}).get("delete_trimming", True))
        ).lower(),
        sortmerna_fq1 = os.path.join("{outdir}", "sortmerna", "{sample_id}.nonrna_1.fastq.gz"),
        sortmerna_fq2 = os.path.join("{outdir}", "sortmerna", "{sample_id}.nonrna_2.fastq.gz"),
        delete_sortmerna = lambda wildcards: str(
            is_enabled("sortmerna") and
            as_bool(config.get("sortmerna", {}).get("delete_sortmerna", True))
        ).lower()
    threads: 1
    resources:
        mem_mb = 1024
    message:
        "{wildcards.sample_id}: Cleaning up temporary files"
    log:
        os.path.join("{outdir}", "logs", "cleanup", "{sample_id}.cleanup.log")
    shell:
        """
        mkdir -p $(dirname {log})
        echo "[INFO] Starting cleanup for {wildcards.sample_id}" > {log}

        # Remove trimmed FASTQ files (conditional)
        if [ "{params.delete_trimming}" = "true" ]; then
            for f in "{params.fq1}" "{params.fq2}" "{params.fq1_fail}" "{params.fq2_fail}"; do
                if [ -f "$f" ]; then
                    rm -f "$f" && echo "[INFO] Removed $f" >> {log}
                fi
            done
        else
            echo "[INFO] delete_trimming=false, skipping trimmed FASTQ deletion." >> {log}
        fi

        # Remove SortMeRNA non-rRNA FASTQ files (conditional)
        if [ "{params.delete_sortmerna}" = "true" ]; then
            for f in "{params.sortmerna_fq1}" "{params.sortmerna_fq2}"; do
                if [ -f "$f" ]; then
                    rm -f "$f" && echo "[INFO] Removed $f" >> {log}
                fi
            done
        else
            echo "[INFO] delete_sortmerna=false, skipping non-rRNA FASTQ deletion." >> {log}
        fi

        # Remove merged raw FASTQs
        for f in "{params.raw_fq1}" "{params.raw_fq2}"; do
            if [ -e "$f" ] || [ -L "$f" ]; then
                rm -f "$f" && echo "[INFO] Removed $f" >> {log}
            fi
        done

        # Remove raw_merged directory if empty
        if [ -d "{params.raw_dir}" ]; then
            if rmdir "{params.raw_dir}" 2>/dev/null; then
                echo "[INFO] Removed empty directory: {params.raw_dir}" >> {log}
            else
                echo "[INFO] Directory not empty, keeping: {params.raw_dir}" >> {log}
            fi
        fi

        # Remove STAR temporary directories and files
        rm -rf "{params.bam_dir}/{wildcards.sample_id}._STARgenome" \
            "{params.bam_dir}/{wildcards.sample_id}._STARpass1" \
            "{params.bam_dir}/{wildcards.sample_id}._STARtmp" 2>/dev/null || true

        find "{params.bam_dir}" \( -name "Aligned.out.sam" \
                -o -name "*.Log.progress.out*" \
                -o -name "*.Log.out*" \
                -o -name "*.SJ.out.tab.bak" \) -print -delete 2>/dev/null || true
        echo "[INFO] Removed STAR temporary files." >> {log}

        # Remove Arriba STAR intermediate files
        find "{wildcards.outdir}/arriba" -name "*.fusionmapping.*" -delete 2>/dev/null || true
        echo "[INFO] Removed Arriba intermediate files." >> {log}

        echo "[INFO] Cleanup completed for {wildcards.sample_id}" >> {log}
        echo "[INFO] Deletion completed for {wildcards.sample_id}. See {log}" > {output.log}
        """
