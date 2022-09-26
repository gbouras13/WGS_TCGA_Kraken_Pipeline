"""
Snakefile for downloading the prebuilt db. Only need to run this once.

snakemake -c 1 -s DownloadDB.smk

"""
import os

# load default config
configfile: os.path.join(workflow.basedir,  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

DBDIR = "Databases"

if not os.path.exists(os.path.join(DBDIR)):
    os.makedirs(os.path.join(DBDIR))

if not os.path.exists(os.path.join(DBDIR, 'standard')):
    os.makedirs(os.path.join(DBDIR, 'standard'))

rule all:
    input:
        os.path.join(DBDIR, 'standard', 'download.dlflag'),
        os.path.join(DBDIR, 'standard' , 'untar.dlflag')

# use the prebuilt standard db
# # https://benlangmead.github.io/aws-indexes/k2

rule download_taxonomy_and_taxs:
    """Rule to Download Pre-built standard db."""
    output:
        os.path.join(DBDIR, 'standard', 'download.dlflag'),
        os.path.join(DBDIR, 'standard', 'k2_pluspf_20210517.tar.gz')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        wget -c "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20210517.tar.gz" -O k2_pluspf_20210517.tar.gz
        mv "k2_pluspf_20210517.tar.gz" {output[1]}
        touch {output[0]}
        """

rule untar:
    """Untar DB File."""
    input:
        os.path.join(DBDIR, 'standard','k2_pluspf_20210517.tar.gz'),
        os.path.join(workflow.basedir, DBDIR,  'standard')
    output:
        os.path.join(DBDIR, 'standard', 'untar.dlflag')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        tar -xzf {input[0]} -C {input[1]}
        touch {output[0]}
        """

