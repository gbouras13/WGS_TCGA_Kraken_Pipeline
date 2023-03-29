"""
Consistent output directory locations
"""


DBDIR = 'Databases'
KRAKENTOOLSDIR = 'Kraken_Tools'
CONTAMINANTS = 'contaminants'

### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'wgs_tcga_out'
else:
    OUTPUT = config['Output']


### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')
WORKDIR = os.path.join(OUTPUT, 'PROCESSING')

# fastq dirs
UNALIGNED_FASTQ  = os.path.join(WORKDIR, 'UNALIGNED_FASTQ')
HOST_UNMAPPED_FASTQ  = os.path.join(WORKDIR, 'HOST_UNMAPPED_FASTQ')
HOST_UNMAPPED_FASTQ_FASTP = os.path.join(WORKDIR, 'HOST_UNMAPPED_FASTQ_FASTP')
BACT_FASTQ = os.path.join(WORKDIR, 'BACT_FASTQ')
VIR_FASTQ = os.path.join(WORKDIR, 'VIR_FASTQ')
FUNGI_FASTQ = os.path.join(WORKDIR, 'FUNGI_FASTQ')
CONCAT_FASTQ = os.path.join(WORKDIR, 'CONCAT_FASTQ')


# dir for flags
LOGS = os.path.join(OUTPUT, 'LOGS')
BIOM = os.path.join(RESULTS, 'BIOM')

# kraken dirs 
KRAKEN = os.path.join(RESULTS, 'KRAKEN')

# bracken dirs 
BRACKEN = os.path.join(RESULTS, 'BRACKEN') 


# get readcount of bams
READCOUNT = os.path.join(RESULTS, 'READCOUNT')





