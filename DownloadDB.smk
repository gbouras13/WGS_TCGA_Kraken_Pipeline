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


### prebuild db
### https://lomanlab.github.io/mockcommunity/mc_databases.html
#
#
#
#

#
# rule all:
#     input:
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree', 'hash.dlflag'),
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree' , 'opts.dlflag'),
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree' , 'taxo.dlflag')
#
#
#
# ### prebuild db
# ### https://lomanlab.github.io/mockcommunity/mc_databases.html
#
#
# rule download_taxonomy_and_taxs:
#     """Generic rule to download a DB file."""
#     output:
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree', 'hash.dlflag'),
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree', 'hash.k2d'),
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree', 'opts.dlflag'),
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree', 'opts.k2d'),
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree', 'taxo.dlflag'),
#         os.path.join(DBDIR, 'kraken2-microbial-fatfree', 'taxo.k2d')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         wget -c "https://refdb.s3.climb.ac.uk/kraken2-microbial/hash.k2d" -O hash.k2d
#         mv "hash.k2d" {output[1]}
#         touch {output[0]}
#         wget "https://refdb.s3.climb.ac.uk/kraken2-microbial/opts.k2d" -O opts.k2d
#         mv "hash.k2d" {output[3]}
#         touch {output[2]}
#         wget "https://refdb.s3.climb.ac.uk/kraken2-microbial/taxo.k2d" -O taxo.k2d
#         mv "hash.k2d" {output[5]}
#         touch {output[4]}
#         """


#############################################
########## old db
#############################################

# # targets
# allDbFiles = []
# for f in config['dbFiles']:
#     allDbFiles.append(os.path.join(DBDIR, f))
#
#
#
# if not os.path.exists(os.path.join(DBDIR, 'bacteria')):
#     os.makedirs(os.path.join(DBDIR, 'bacteria'))
#
# if not os.path.exists(os.path.join(DBDIR, 'bacteria', 'taxonomy')):
#     os.makedirs(os.path.join(DBDIR, 'bacteria', 'taxonomy'))
#
#
# rule all:
#     input:
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'taxdump.dlflag'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'taxdump.tarflag'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'nucl_gb.dlflag'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'nucl_wgs.dlflag'),
#         os.path.join(DBDIR, 'bacteria', 'build.dlflag')
#
# rule download_taxonomy_and_taxs:
#     """Generic rule to download a DB file."""
#     output:
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'taxdump.dlflag'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'taxdump.tar.gz'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'nucl_gb.dlflag'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'nucl_gb.accession2taxid.gz'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'nucl_wgs.dlflag'),
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'nucl_wgs.accession2taxid.gz')
#     conda:
#         os.path.join('envs','kraken2.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         wget "ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" -O taxdump.tar.gz
#         mv "taxdump.tar.gz" {output[2]}
#         touch {output[1]}
#         wget "ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz" -O nucl_gb.accession2taxid.gz
#         mv "nucl_gb.accession2taxid.gz" {output[4]}
#         touch {output[3]}
#         wget "ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz" -O nucl_wgs.accession2taxid.gz
#         mv "nucl_wgs.accession2taxid.gz" {output[6]}
#         touch {output[5]}
#         """
#
# rule untar:
#     """Generic rule to download a DB file."""
#     input:
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'taxdump.tar.gz'),
#         os.path.join(workflow.basedir, DBDIR,  'bacteria', 'taxonomy')
#     output:
#         os.path.join(DBDIR, 'bacteria', 'taxonomy' , 'taxdump.tarflag')
#     conda:
#         os.path.join('envs','kraken2.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         tar zxf {input[0]} -C {input[1]}
#         touch {output[0]}
#         # kraken2-build --download-library bacteria --db {output[0]}
#         """
#
# ####### https://github.com/DerrickWood/kraken2/issues/518
#
# rule build:
#     """build."""
#     input:
#         os.path.join(workflow.basedir, DBDIR, 'bacteria')
#     output:
#         os.path.join(DBDIR, 'bacteria', 'build.dlflag')
#     conda:
#         os.path.join('envs','kraken2.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         kraken2-build --download-library bacteria --db {input[0]}
#         touch {output[0]}
#         """
#
#
#
#
# # rule build_db_file:
# #     """Generic rule to download a DB file."""
# #     output:
# #         os.path.join(DBDIR, 'bacteria')
# #     conda:
# #         os.path.join('envs','kraken2.yaml')
# #     threads:
# #         BigJobCpu
# #     resources:
# #         mem_mb=BigJobMem
# #     shell:
# #         """
# #         kraken2-build --build --db {output[0]} --threads 8
# #         """
