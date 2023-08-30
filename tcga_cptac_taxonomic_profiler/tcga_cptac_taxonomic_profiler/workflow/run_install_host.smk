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

rule get_host:
    """ 
    This can definitely be improved 
    """
    params:
        host_db = HostDir
    conda:
        os.path.join( 'envs', 'gzip.yaml')
    output:
        fasta = os.path.join(HostDir,"human-t2t-hla.fa")
    shell:
        """
        cd {params.host_db}
        wget "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz" -O human-t2t-hla.fa.gz
        gunzip human-t2t-hla.fa.gz
        """

rule get_phix:
    """ 
    This can definitely be improved 
    """
    params:
        host_db = HostDir
    conda:
        os.path.join( 'envs', 'gzip.yaml')
    output:
        fasta = os.path.join(HostDir,"human-t2t-hla.fa"),
        fasta_phix = os.path.join(HostDir,"phix.fna")
    shell:
        """
        cd {params.host_db}
        wget ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Sinsheimervirus_phiX174/latest_assembly_versions/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz   -O phix.fna.gz
        gunzip phix.fna.gz 
        """

rule combine_phix174:
    """ 
    
    """
    input:
        chm13 = os.path.join(HostDir,"human-t2t-hla.fa"),
        phix174 = os.path.join(HostDir,"phix.fna")
    output:
        fasta = os.path.join(HostDir,"human-t2t-hla-phix174.fa")
    shell:
        """
        cat {input.chm13} {input.phix174} > {output.fasta}
        """