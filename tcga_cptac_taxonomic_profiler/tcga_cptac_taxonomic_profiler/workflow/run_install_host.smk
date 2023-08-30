"""
Snakefile for downloading CHM
"""
import os
import attrmap as ap
import attrmap.utils as au

# load default config
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
config = ap.AttrMap(config)



# config

HostDir = config.databases.host_db

if not os.path.exists(os.path.join(HostDir)):
    os.makedirs(os.path.join(HostDir))

rule all:
    input:
        os.path.join(HostDir,"human-t2t-hla-phix174.fa")

rule get_db:
    """ 
    This can definitely be improved 
    """
    params:
        host_db = HostDir
    conda:
        os.path.join( 'envs', 'gzip.yml')
    output:
        fasta = os.path.join(HostDir,"human-t2t-hla.fa"),
        fasta_phix = os.path.join(HostDir,"NC_001422.fna")
    shell:
        """
        cd {params.host_db}
        wget "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz" -O human-t2t-hla.fa.gz
        gunzip human-t2t-hla.fa.gz
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna  -O NC_001422.fna
        """

rule combine_phix174:
    """ 
    
    """
    input:
        chm13 = os.path.join(HostDir,"human-t2t-hla.fa"),
        phix174 = os.path.join(HostDir,"NC_001422.fna")
    output:
        fasta = os.path.join(HostDir,"human-t2t-hla-phix174.fa")
    shell:
        """
        cat {input.chm13} {input.phix174} > {output.fasta}
        """