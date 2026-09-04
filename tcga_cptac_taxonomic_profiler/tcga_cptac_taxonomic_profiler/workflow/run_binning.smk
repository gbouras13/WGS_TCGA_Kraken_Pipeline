"""
The snakefile that runs the binning pipeline.

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
include: "rules/sample_assembly_binning.smk"
include: "rules/multi_binners.smk"

# import targets
include: "rules/targets.smk"

rule all:
    input:
        SampleAssemblyBins

