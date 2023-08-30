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
        fasta = os.path.join(HostDir,"human-t2t-hla.fa")
    shell:
        """
        cd {params.host_db}
        wget "https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz" -O human-t2t-hla.fa.gz
        gunzip human-t2t-hla.fa.gz
        """

rule combine_phix174:
    """ 
    
    """
    input:
        chm13 = os.path.join(HostDir,"human-t2t-hla.fa"),
        phix174 = os.path.join( 'db',"phix174.fasta")
    output:
        fasta = os.path.join(HostDir,"human-t2t-hla-phix174.fa")
    shell:
        """
        cat {input.chm13} {input.phix174} > {output.fasta}
        """