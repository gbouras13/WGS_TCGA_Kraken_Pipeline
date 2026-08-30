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
HOST_RM_FASTQ    = os.path.join(RESULTS, 'HOST_RM_FASTQ')

# dir for flags
FLAGS = os.path.join(OUTPUT, 'FLAGS')
BIOM = os.path.join(RESULTS, 'BIOM')

# kraken and bracken dirs 
KRAKEN = os.path.join(RESULTS, 'KRAKEN')
BRACKEN = os.path.join(RESULTS, 'BRACKEN') 

# metabuli (joint AA+DNA classifier, independent of kraken2)
METABULI = os.path.join(RESULTS, 'METABULI')

# sylph (ANI-based containment profiler)
SYLPH = os.path.join(RESULTS, 'SYLPH')

# mmseqs2 
FASTA = os.path.join(PROCESSING, 'FASTA')
MMSEQS2 = os.path.join(RESULTS, 'MMSEQS2')

# get readcount of bams
READCOUNT = os.path.join(RESULTS, 'READCOUNT')

# assembly
SAMPLE_ASSEMBLIES = os.path.join(RESULTS, 'SAMPLE_ASSEMBLIES')

# singlem 
SINGLEM = os.path.join(RESULTS, 'SINGLEM')

#binning 
# sample
VAMB_CATALOGUE = os.path.join(RESULTS, 'VAMB_CATALOGUE')
VAMB_BAMS = os.path.join(RESULTS, 'VAMB_BAMS')
VAMB_RESULTS = os.path.join(RESULTS, 'VAMB_RESULTS')
SEMIBIN2_RESULTS =   os.path.join(RESULTS, 'SEMIBIN2_RESULTS')
METABAT2_RESULTS =   os.path.join(RESULTS, 'METABAT2_RESULTS')
CONCOCT_RESULTS  =   os.path.join(RESULTS, 'CONCOCT_RESULTS')
COMEBIN_RESULTS  =   os.path.join(RESULTS, 'COMEBIN_RESULTS')
BINETTE_RESULTS  =   os.path.join(RESULTS, 'BINETTE_RESULTS')
ALL_MAGS = os.path.join(RESULTS, 'ALL_MAGS')

# checkm2
CHECKM2_RESULTS  = os.path.join(RESULTS, 'CHECKM2_RESULTS')

# GTDB
GTDB_OUTDIR = os.path.join(RESULTS, 'GTDB_RESULTS')
GTDB_MASH_OUTDIR = os.path.join(GTDB_OUTDIR, 'MASH')

# BAKTA
BAKTA = os.path.join(RESULTS, 'BAKTA')

# PHISPY
PHISPY = os.path.join(RESULTS, 'PHISPY')

# MGE
GENOMAD_RESULTS = os.path.join(RESULTS, 'GENOMAD_RESULTS')