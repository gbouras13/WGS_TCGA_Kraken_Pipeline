"""
The snakefile that extracts the unaligned fastqs files frome the raw TCGA/CPTAC bams

Also extracts the read counts for the bams in total


"""

import glob
import attrmap as ap
import attrmap.utils as au


configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
config = ap.AttrMap(config)

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
BigJobMem = config.BigJobMem
BigJobCpu = config.BigJobCpu
MediumJobMem = config.MediumJobMem
SmallJobMem = config.SmallJobMem


### DIRECTORIES
# get if needed
BAMS_DIR = config.input
OUTPUT = config.output
THREADS = config.threads

### DIRECTORIES
include: "rules/directories.smk"

# get input and put directorr

# Parse the samples from *.bam files

include: "rules/samples_from_bam.smk"
sampleReads = parseSamplesBam(BAMS_DIR)
SAMPLES = sampleReads.keys()

# Import rules and functions
include: "rules/targets.smk"
include: "rules/bam_to_fastq.smk"
include: "rules/read_counts.smk"

rule all:
    input:
        ExtractFiles
