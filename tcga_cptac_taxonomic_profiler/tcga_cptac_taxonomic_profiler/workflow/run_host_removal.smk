"""
Deplete human reads from the extracted unaligned FASTQs with Deacon.

Runs between `extract` and `kraken`/`singlem`/`assembly`. Previously this step
was performed outside the workflow, so the parameters used were not recorded.
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
config = ap.AttrMap(config)

### DIRECTORIES
INPUT = config.input
OUTPUT = config.output
THREADS = config.threads

include: "rules/directories.smk"

# Parse samples from the extracted unaligned FASTQs
SAMPLES, = glob_wildcards(os.path.join(INPUT, "{sample}_R1.fastq.gz"))

include: "rules/targets.smk"
include: "rules/host_removal.smk"

rule all:
    input:
        HostRemovalTargets
