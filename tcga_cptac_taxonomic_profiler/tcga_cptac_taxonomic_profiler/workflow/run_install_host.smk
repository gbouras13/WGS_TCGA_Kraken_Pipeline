"""
Snakefile for downloading CHM
"""
import os
import attrmap as ap
import attrmap.utils as au

# load default config
# Ignore ~/.local/lib/pythonX.Y/site-packages in every rule.
#
# A user-site numpy 1.26.4 under /home shadowed the conda env's numpy 2.5.2 in
# COMEBin, whose scipy 1.18 requires numpy >= 2. The env was built correctly;
# python simply preferred the home directory. It fails as
# "module 'numpy' has no attribute 'long'", which reads like a version pin
# problem in the env and is not - nothing in the env is wrong.
#
# This is set here rather than per-rule because it can affect ANY environment
# whose interpreter version matches something installed under ~/.local.
shell.prefix("export PYTHONNOUSERSITE=1; ")

configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
config = ap.AttrMap(config)



# config

HostDir = config.databases.host_db

if not os.path.exists(os.path.join(HostDir)):
    os.makedirs(os.path.join(HostDir))

rule all:
    input:
        os.path.join(HostDir,"human-t2t-hla-phix174.fa"),
        os.path.join(HostDir,"panhuman-1.k31w15.idx")

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

rule get_deacon_index:
    """Fetch the prebuilt Deacon panhuman-1 index (~3.3 GB).

    panhuman-1 = HPRC year-1 assemblies + CHM13v2.0 + GRCh38.p14, with
    bacterial and viral sequence removed. `deacon index fetch` resolves the
    current mirror itself rather than hardcoding a URL.
    Zenodo archive for the record: https://zenodo.org/records/17288185
    """
    params:
        host_db = HostDir
    conda:
        os.path.join('envs', 'deacon.yaml')
    output:
        index = os.path.join(HostDir, "panhuman-1.k31w15.idx")
    shell:
        """
        deacon index fetch panhuman-1 --output {output.index}
        """
