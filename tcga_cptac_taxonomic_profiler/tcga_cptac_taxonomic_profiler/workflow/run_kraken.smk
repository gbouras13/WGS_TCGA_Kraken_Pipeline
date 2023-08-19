"""
The snakefile that runs the kraken pipeline.

tcga_cptac_taxonomic_profiler extract needs to be run first.

"""


import glob


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



### DIRECTORIES
include: "rules/directories.smk"

# get if needed
INPUT = config['input']
OUTPUT = config['output']
KRAKENDB = config["kraken_db"]

# Parse the samples with metasnek

# https://gist.github.com/beardymcjohnface/bb161ba04ae1042299f48a4849e917c8

# this parses the samples into a dictionary
sample_dict = fastq_finder.parse_samples_to_dictionary(INPUT)
SAMPLES = list(sample_dict.keys())

# Import rules and functions
include: "rules/targets.smk"
include: "rules/kraken.smk"
include: "rules/bracken.smk"
include: "rules/biom.smk"



rule all:
    input:
        KrakenTargets

