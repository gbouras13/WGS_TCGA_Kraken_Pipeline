"""
The snakefile that extracts the unaligned fastqs files frome the raw TCGA bams

Also extracts the read counts for the bams in total

Manual launch example:

snakemake -c 1 -s extract_unaligned_fastq.smk --use-conda --config Bams=Bams/  --conda-create-envs-only --conda-frontend conda

compute node

snakemake -c 16 -s extract_unaligned_fastq.smk --use-conda --config Bams=Bams/ Output=out/

requires Bams as input - the file containing all the TCGA Bams

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
BAMS = config['Bams']
OUTPUT = config['Output']

# Parse the samples and read files
# requires 'Reads' as sample input
include: "rules/samples_from_bam.smk"
sampleReads = parseSamplesBam(BAMS)
SAMPLES = sampleReads.keys()

# Import rules and functions
include: "rules/targets.smk"
include: "rules/bam_to_fastq.smk"
include: "rules/read_counts.smk"

rule all:
    input:
        BamToFastqFiles


