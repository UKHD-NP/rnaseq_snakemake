from urllib.request import urlretrieve
import gzip
import shutil
import ssl
import os
import sys
ssl._create_default_https_context = ssl._create_unverified_context

def info(msg):
    print(f"[INFO] {msg}")


def warning(msg):
    print(f"[WARNING] {msg}")


def fatal(msg):
    sys.exit(f"[ERROR] {msg}")


REFERENCE_URLS = {
    "hg19": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz",
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
    ),
    "hg38": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz",
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.basic.annotation.gtf.gz",
    ),
    "chm13v2": (
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
        "https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz",
    ),
    "m39": (
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.primary_assembly.genome.fa.gz",
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.primary_assembly.annotation.gtf.gz",
    ),
}

assembly = config["ref"]["assembly"]
custom_path_fasta = None
custom_path_gtf = None
url_fasta = None
url_gtf = None

if assembly == "custom":
    custom_path_fasta = config["ref"].get("custom_fasta_path")
    custom_path_gtf = config["ref"].get("custom_gtf_path")
    if not custom_path_fasta or not custom_path_gtf:
        fatal("For custom assembly, set both 'ref.custom_fasta_path' and 'ref.custom_gtf_path' in config/config.yaml.")
elif assembly in REFERENCE_URLS:
    url_fasta, url_gtf = REFERENCE_URLS[assembly]
else:
    fatal(f"Invalid 'ref.assembly': '{assembly}'. Allowed values: hg19, hg38, chm13v2, m39, custom.")

# Set paths to reference files based on organism and assembly selection
ref_dir = os.path.join("references", config["ref"]["assembly"])
os.makedirs(ref_dir, exist_ok=True)
config["ref"]["fasta"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".fa")
config["ref"]["gtf"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".gtf")
config["ref"]["tx_fasta"] = os.path.join(ref_dir, config["ref"]["assembly"] + "_tx.fa")
config["ref"]["star_idx"] = os.path.join(ref_dir, config["ref"]["assembly"] + "_staridx")

# Set BED file path if not provided or empty (for RSeQC modules)
if not config["ref"].get("bed") or config["ref"]["bed"] == "":
    config["ref"]["bed"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".bed")

# Override STAR index path if custom path is provided
if "custom_star_index_path" in config["ref"] and config["ref"]["custom_star_index_path"]:
    custom_star_idx = config["ref"]["custom_star_index_path"]
    if os.path.isdir(custom_star_idx):
        info(f"Using custom STAR index from: {custom_star_idx}")
        config["ref"]["star_idx"] = custom_star_idx
    else:
        warning(f"Custom STAR index path does not exist: {custom_star_idx}. STAR index will be generated.")

# Helper function to decompress files
def decompress_file(input_file, output_file):
    """Decompress gzipped file to output path with error handling"""
    try:
        with gzip.open(input_file, "rb") as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return True
    except Exception as e:
        warning(f"Failed to decompress '{input_file}': {str(e)}")
        return False

def stage_local_reference(source_path, target_path, label):
    """Copy or decompress a local reference file into the references directory."""
    if source_path.endswith(".gz"):
        info(f"Decompressing custom {label} file: {source_path}")
        if not decompress_file(source_path, target_path):
            fatal(f"Could not decompress custom {label} file: {source_path}")
    else:
        try:
            shutil.copyfile(source_path, target_path)
        except Exception as e:
            fatal(f"Could not copy custom {label} file '{source_path}' -> '{target_path}': {str(e)}")

def download_reference_file(url, compressed_path, target_path, label):
    """Download and decompress a remote reference file."""
    try:
        info(f"Downloading {label} reference from: {url}")
        urlretrieve(url, compressed_path)
        info(f"Decompressing downloaded {label} file: {compressed_path}")
        if not decompress_file(compressed_path, target_path):
            fatal(f"Could not decompress downloaded {label} file: {compressed_path}")
        os.remove(compressed_path)
    except Exception as e:
        fatal(f"Could not download {label} file from '{url}': {str(e)}")

# Check if FASTA reference already exists. Download it if not.
if not os.path.isfile(config["ref"]["fasta"]):
    if config["ref"]["assembly"] == "custom":
        info("Using custom FASTA reference.")
        stage_local_reference(custom_path_fasta, config["ref"]["fasta"], "FASTA")
    else:
        info(f"Downloading {config['ref']['assembly']} FASTA reference.")
        fasta_compressed = os.path.join(ref_dir, config["ref"]["assembly"] + ".fa.gz")
        download_reference_file(url_fasta, fasta_compressed, config["ref"]["fasta"], "FASTA")

# Check if GTF annotation already exists. Download it if not.
if not os.path.isfile(config["ref"]["gtf"]):
    if config["ref"]["assembly"] == "custom":
        info("Using custom GTF annotation.")
        stage_local_reference(custom_path_gtf, config["ref"]["gtf"], "GTF")
    else:
        info(f"Downloading {config['ref']['assembly']} GTF annotation.")
        gtf_compressed = os.path.join(ref_dir, config["ref"]["assembly"] + ".gtf.gz")
        download_reference_file(url_gtf, gtf_compressed, config["ref"]["gtf"], "GTF")

rule samtools_faidx:
    # Generate FASTA index
    input:
        config['ref']['fasta']
    output:
        config['ref']['fasta'] + ".fai"
    message:
        "Generating FASTA index"
    conda:
        os.path.join(workflow.basedir, "envs", "samtools.yml")
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
    message:
        "Generating transcript fasta"
    params:
        is_custom = "true" if (config['ref']['assembly'] == "custom" and 'custom_transcript_path' in config['ref']) else "false",
        custom_path = config['ref'].get('custom_transcript_path', ""),
        is_gz = "true" if (config['ref'].get('custom_transcript_path', "").endswith('.gz')) else "false"
    conda:
        os.path.join(workflow.basedir, "envs", "gffread.yml")
    threads: 1
    log:
        os.path.join(ref_dir, "gffread", f"{config['ref']['assembly']}.log")
    shell:
        """
        if [ "{params.is_custom}" = "true" ] && [ -f "{params.custom_path}" ]; then
            echo "[INFO] Using provided custom transcript file: {params.custom_path}"
            if [ "{params.is_gz}" = "true" ]; then
                echo "[INFO] Decompressing custom transcript file..."
                gunzip -c "{params.custom_path}" > "{output}" || {{ echo "[ERROR] Failed to decompress custom transcript file: {params.custom_path}"; exit 1; }}
            else
                echo "[INFO] Copying custom transcript file..."
                cp "{params.custom_path}" "{output}" || {{ echo "[ERROR] Failed to copy custom transcript file: {params.custom_path}"; exit 1; }}
            fi
        else
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
        fi
        """

rule gtf2bed:
    # Generate BED file from GTF for RSeQC modules
    input:
        config['ref']['gtf']
    output:
        config['ref']['bed']
    message:
        "Converting GTF to BED format for RSeQC"
    params:
        gtf2bed_script = os.path.join(workflow.basedir, "scripts", "gtf2bed")
    threads: 1
    log:
        os.path.join(ref_dir, "gtf2bed", f"{config['ref']['assembly']}.log")
    shell:
        """
        # Use gtf2bed script if available, otherwise use awk
        if [ -f "{params.gtf2bed_script}" ]; then
            echo "[INFO] Using gtf2bed script: {params.gtf2bed_script}"
            "{params.gtf2bed_script}" {input} > {output} 2> {log}
        else
            echo "[WARNING] gtf2bed script not found, using awk fallback."
            awk 'BEGIN{{OFS="\\t"}} $3=="exon" {{print $1,$4-1,$5,$10,$6,$7,$4-1,$5,"0,0,0",1,$5-$4+1,0}}' {input} > {output} 2> {log}
        fi

        # Check if output was created successfully
        if [ ! -s "{output}" ]; then
            echo "[ERROR] Failed to create BED file: {output}"
            exit 1
        fi

        echo "[INFO] BED file created successfully: {output}"
        """
