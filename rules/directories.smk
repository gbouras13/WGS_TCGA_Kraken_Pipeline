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
BACT_FASTQ_FIRST_PASS = os.path.join(WORKDIR, 'BACT_FASTQ_FIRST_PASS')
VIR_FASTQ_FIRST_PASS = os.path.join(WORKDIR, 'VIR_FASTQ_FIRST_PASS')
FUNGI_FASTQ_FIRST_PASS = os.path.join(WORKDIR, 'FUNGI_FASTQ_FIRST_PASS')
CONCAT_FASTQ = os.path.join(WORKDIR, 'CONCAT_FASTQ')

# dir for flags
LOGS = os.path.join(OUTPUT, 'LOGS')
BIOM = os.path.join(RESULTS, 'BIOM')

# kraken dirs 

KRAKEN_FIRST_PASS = os.path.join(RESULTS, 'KRAKEN_FIRST_PASS')
KRAKEN_SECOND_PASS = os.path.join(RESULTS, 'KRAKEN_SECOND_PASS')

# bracken dirs 

BRACKEN_FIRST_PASS = os.path.join(RESULTS, 'BRACKEN_FIRST_PASS') 
BRACKEN_SECOND_PASS = os.path.join(RESULTS, 'BRACKEN_SECOND_PASS') 

# MEGAHIT
MEGAHIT = os.path.join(RESULTS, 'MEGAHIT')
# needs to be created before megahit is run for some reason
if not os.path.exists(MEGAHIT):
  os.makedirs(MEGAHIT)



# get readcount of bams
READCOUNT = os.path.join(RESULTS, 'READCOUNT')





