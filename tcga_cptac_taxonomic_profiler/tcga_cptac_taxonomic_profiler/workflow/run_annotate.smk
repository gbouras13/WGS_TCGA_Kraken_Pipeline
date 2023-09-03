"""
The snakefile that runs the binning pipeline.

"""

import glob
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
# get all the mags
include: "rules/mags_from_magdir.smk"
mag_dict = parseSamplesMags(ALL_MAGS)
HQ_MED_MAGS = mag_dict.keys()


# Import rules and functions
include: "rules/bakta.smk"
include: "rules/phispy.smk"

# import targets
include: "rules/targets.smk"

rule all:
    input:
        AnnotationTargets

