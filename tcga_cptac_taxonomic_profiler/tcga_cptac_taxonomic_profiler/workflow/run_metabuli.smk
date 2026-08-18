"""
Classify host-depleted reads with Metabuli (joint amino acid + DNA).

Runs on the same input as the kraken stage, as an algorithmically independent
second profiler. Compare the two rather than choosing between them.
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
config = ap.AttrMap(config)

### DIRECTORIES
INPUT = config.input
OUTPUT = config.output
THREADS = config.threads

include: "rules/directories.smk"

SAMPLES, = glob_wildcards(os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz"))

include: "rules/targets.smk"
include: "rules/metabuli.smk"

rule all:
    input:
        MetabuliTargets
