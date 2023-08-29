"""
Consistent output directory locations
"""


BENCHMARKS = os.path.join(OUTPUT, 'BENCHMARKS')
LOGS = os.path.join(OUTPUT, 'LOGS')

### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')
PROCESSING = os.path.join(OUTPUT, 'PROCESSING')

# fastq dirs
UNALIGNED_FASTQ  = os.path.join(RESULTS, 'UNALIGNED_FASTQ')

# dir for flags
FLAGS = os.path.join(OUTPUT, 'FLAGS')
BIOM = os.path.join(RESULTS, 'BIOM')

# kraken and bracken dirs 
KRAKEN = os.path.join(RESULTS, 'KRAKEN')
BRACKEN = os.path.join(RESULTS, 'BRACKEN') 

# mmseqs2 
FASTA = os.path.join(PROCESSING, 'FASTA')
MMSEQS2 = os.path.join(RESULTS, 'MMSEQS2')

# get readcount of bams
READCOUNT = os.path.join(RESULTS, 'READCOUNT')

# assembly
SAMPLE_ASSEMBLIES = os.path.join(RESULTS, 'SAMPLE_ASSEMBLIES')
COASSEMBLY = os.path.join(RESULTS, 'COASSEMBLY')
COASSEMBLY_RESULTS = os.path.join(RESULTS, 'COASSEMBLY_RESULTS')

#binning 
# sample
VAMB_CATALOGUE = os.path.join(RESULTS, 'VAMB_CATALOGUE')
VAMB_BAMS = os.path.join(RESULTS, 'VAMB_BAMS')
VAMB_RESULTS = os.path.join(RESULTS, 'VAMB_RESULTS')
SEMIBIN2_RESULTS =   os.path.join(RESULTS, 'SEMIBIN2_RESULTS')
ALL_MAGS = os.path.join(RESULTS, 'ALL_MAGS')

# checkm2
CHECKM2_RESULTS  = os.path.join(RESULTS, 'CHECKM2_RESULTS')

# GTDB
GTDB_OUTDIR = os.path.join(RESULTS, 'GTDB_RESULTS')
GTDB_MASH_OUTDIR = os.path.join(GTDB_OUTDIR, 'MASH')

# MGE
GENOMAD_RESULTS = os.path.join(RESULTS, 'GENOMAD_RESULTS')