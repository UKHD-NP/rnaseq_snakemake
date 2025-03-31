from urllib.request import urlretrieve
import gzip
import shutil
import subprocess
import ssl
ssl._create_default_https_context = ssl._create_unverified_context


# Set URLs for FASTA and GTF download
url = ""
if config['ref']['assembly'] == "hg19":
    url_fasta = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz"
    url_gtf = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
elif config['ref']['assembly'] == "hg38":
    url_fasta = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz"
    url_gtf = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.basic.annotation.gtf.gz"
elif config['ref']['assembly'] == "chm13v2":
    url_fasta = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    url_gtf = "https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz"
elif config['ref']['assembly'] == "m39":
    url_fasta = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.primary_assembly.genome.fa.gz"
    url_gtf = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.primary_assembly.annotation.gtf.gz"
else:
    sys.exit("Wrong genome assembly name. Exiting.")


# Set paths to reference files based on organism and assembly selection
ref_dir = os.path.join("references", config['ref']['assembly'])
os.makedirs(ref_dir, exist_ok=True)
config['ref']['fasta'] = os.path.join(ref_dir, config['ref']['assembly'] + ".fa")
config['ref']['gtf'] = os.path.join(ref_dir, config['ref']['assembly'] + ".gtf")
config['ref']['tx_fasta'] = os.path.join(ref_dir, config['ref']['assembly'] + "_tx.fa")
config['ref']['star_idx'] = os.path.join(ref_dir, config['ref']['assembly'] + "_staridx")


# Check if FASTA reference already exists. Download it if not.
if not os.path.isfile(config['ref']['fasta']):
    print("Downloading "+config['ref']['assembly']+" FASTA reference...")
    fasta_compressed = os.path.join(ref_dir, config['ref']['assembly'] + ".fa.gz")
    urlretrieve(url_fasta, fasta_compressed)    
    print("Decompressing FASTA file...")
    with gzip.open(fasta_compressed, "rb") as f_in:
        with open(config['ref']['fasta'], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)    
    os.remove(fasta_compressed)


# Check if GTF annotation already exists. Download it if not.
if not os.path.isfile(config['ref']['gtf']):
    print("Downloading "+config['ref']['assembly']+" GTF annotation...")
    gtf_compressed = os.path.join(ref_dir, config['ref']['assembly'] + ".gtf.gz")
    urlretrieve(url_gtf, gtf_compressed)    
    print("Decompressing GTF file...")
    with gzip.open(gtf_compressed, "rb") as f_in:
        with open(config['ref']['gtf'], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)    
    os.remove(gtf_compressed)


rule samtools_faidx:
    input:
        config['ref']['fasta']
    output:
        config['ref']['fasta'] + ".fai"
    message:
        "Generating FASTA index"
    conda:
        "../envs/samtools.yml"
    threads: 1
    shell:
        """
        samtools faidx {input}
        """

rule gffread:
    input:
        fasta = config['ref']['fasta'],
        gtf = config['ref']['gtf']
    output:
        config['ref']['tx_fasta']
    message:
        "Generating transcript fasta"
    conda:
        "../envs/gffread.yml"
    threads: 1
    shell:
        """
        gffread \
            -g {input.fasta} \
            -w {output} \
            {input.gtf}
        """
