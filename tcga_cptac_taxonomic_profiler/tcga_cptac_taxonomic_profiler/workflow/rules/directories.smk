"""
Consistent output directory locations
"""


### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'wgs_tcga_out'
else:
    OUTPUT = config['Output']


BENCHMARKS = os.path.join(OUTPUT, 'BENCHMARKS')
LOGS = os.path.join(OUTPUT, 'LOGS')

### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')

# fastq dirs
UNALIGNED_FASTQ  = os.path.join(RESULTS, 'UNALIGNED_FASTQ')

# dir for flags
FLAGS = os.path.join(OUTPUT, 'FLAGS')
BIOM = os.path.join(RESULTS, 'BIOM')

# kraken and bracken dirs 
KRAKEN = os.path.join(RESULTS, 'KRAKEN')
BRACKEN = os.path.join(RESULTS, 'BRACKEN') 

# get readcount of bams
READCOUNT = os.path.join(RESULTS, 'READCOUNT')





