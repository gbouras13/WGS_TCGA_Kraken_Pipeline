"""
Profile host-depleted reads with sylph (ANI-based containment estimation).

Runs on the same input as the kraken and metabuli stages. Sylph contributes
containment ANI per detected genome, which neither read classifier provides.
"""

import glob
import attrmap as ap
import attrmap.utils as au


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
include: "rules/sylph.smk"

rule all:
    input:
        SylphTargets
