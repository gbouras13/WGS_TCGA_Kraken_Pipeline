"""
The snakefile that runs the pipeline.
Manual launch example:

snakemake -c 1 -s wgs_runner.smk --use-conda --config Reads=Bams/  --conda-create-envs-only --conda-frontend conda
compute node
snakemake -c 16 -s wgs_runner.smk --use-conda --config Reads=Bams/ Output=out/
"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir,  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

### DIRECTORIES
include: "rules/directories.smk"

# get if needed
READS = config['Reads']
OUTPUT = config['Output']

# Parse the samples and read files
include: "rules/samples.smk"
sampleReads = parseSamples(READS)
SAMPLES = sampleReads.keys()

# Import rules and functions
include: "rules/targets.smk"
include: "rules/bam_to_fastq.smk"
include: "rules/kraken.smk"
include: "rules/bracken.smk"
include: "rules/assembly.smk"

rule all:
    input:
        ## Preprocessing
        PreprocessingFiles
        # ## Assembly
        # AssemblyFiles,
        # ## Translated (nt-to-aa) search
        # SecondarySearchFilesAA,
        # ## Untranslated (nt-to-nt) search
        # SecondarySearchFilesNT,
        # ## Contig annotation
        # ContigAnnotFiles,
        # ## Mapping (read-based contig id)
        # MappingFiles,
        # ## Summary
        # SummaryFiles
