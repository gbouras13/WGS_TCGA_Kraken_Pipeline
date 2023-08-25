"""
The snakefile that runs the assembly pipeline.

tcga_cptac_taxonomic_profiler extract needs to be run first.

"""

import glob
from metasnek import fastq_finder
import attrmap as ap
import attrmap.utils as au


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
# from hecatomb
config = ap.AttrMap(config)


### DIRECTORIES
# get if needed
INPUT = config.input
OUTPUT = config.output
TMPDIR = config.tmpdir


include: "rules/directories.smk"

# this parses the samples into a dictionary
sample_dict = fastq_finder.parse_samples_to_dictionary(INPUT)
SAMPLES = list(sample_dict.keys())


# Import rules and functions
include: "rules/co_assembly.smk"
include: "rules/sample_assembly.smk"

# import targets
include: "rules/targets.smk"

rule all:
    input:
        SampleAssemblyBins

