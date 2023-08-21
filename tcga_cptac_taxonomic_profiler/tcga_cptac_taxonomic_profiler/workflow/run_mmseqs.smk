"""
The snakefile that runs the mmseqs2 pipeline.

tcga_cptac_taxonomic_profiler extract needs to be run first.

"""

import glob
from metasnek import fastq_finder

configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')


# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + config['log'])

onsuccess:
    copy_log_file()

onerror:
    copy_log_file()



### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
MediumJobMem = config["MediumJobMem"]
SmallJobMem = config["SmallJobMem"]
BigJobTimeMin = config["BigJobTimeMin"]
MediumJobTimeMin = config["MediumJobTimeMin"]
SmallJobTimeMin = config['SmallJobTimeMin']

### DIRECTORIES
# get if needed
INPUT = config['input']
OUTPUT = config['output']
MMSEQS2_DB = config['mmseqs2_db']

include: "rules/directories.smk"

# this parses the samples into a dictionary
sample_dict = fastq_finder.parse_samples_to_dictionary(INPUT)
SAMPLES = list(sample_dict.keys())


# Import rules and functions
include: "rules/fastq_to_fasta.smk"
include: "rules/mmseqs2_easy_tax.smk"



rule all:
    input:
        MMseqsTargets

