"""
The snakefile that runs the pipeline.

extract_unaligned_fastq.smk needs to be run first.

snakemake -c 1 -s wgs_runner.smk --use-conda --config Reads=Bams/  --conda-create-envs-only --conda-frontend conda
compute node
snakemake -c 16 -s wgs_runner.smk --use-conda --config Output=out/

"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir,  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
MediumJobMem = config["MediumJobMem"]
SmallJobMem = config["SmallJobMem"]

### DIRECTORIES
include: "rules/directories.smk"

# get if needed
OUTPUT = config['Output']

# Parse the samples and read files
include: "rules/samples_from_output.smk"
sampleReads = parseSamples(UNALIGNED_FASTQ)
SAMPLES = sampleReads.keys()

# Import rules and functions
include: "rules/targets.smk"
include: "rules/kraken_first_pass.smk"
include: "rules/extract_bact_vir_reads.smk"
include: "rules/fastp.smk"
include: "rules/kraken_second_pass.smk"
include: "rules/bracken.smk"



rule all:
    input:
        OutputFiles

