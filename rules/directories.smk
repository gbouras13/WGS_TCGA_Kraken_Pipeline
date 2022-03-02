"""
Database and output locations for Hecatomb
Ensures consistent variable names and file locations for the pipeline, the database download script,
and the addHost script.
"""


DBDIR = 'Databases'
KRAKENTOOLSDIR = 'Kraken_Tools'

### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'wgs_tcga_out'
else:
    OUTPUT = config['Output']


### DATABASE SUBDIRs
CONPATH = os.path.join(DBDIR, "contaminants")
TAX = os.path.join(DBDIR, "tax", "taxonomy")
TABLES = os.path.join(DBDIR, "tables")
HOSTPATH = os.path.join(DBDIR, "host")


### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')
WORKDIR = os.path.join(OUTPUT, 'PROCESSING')
TMP = os.path.join(WORKDIR, 'TMP')
LOGS = os.path.join(OUTPUT, 'LOGS')
BIOM = os.path.join(RESULTS, 'BIOM')
KRAKEN_S = os.path.join(RESULTS, 'KRAKEN_S')
KRAKEN_G = os.path.join(RESULTS, 'KRAKEN_G')
BRACKEN = os.path.join(RESULTS, 'BRACKEN')
MEGAHIT = os.path.join(RESULTS, 'MEGAHIT')


# needs to be created before megahit is run
if not os.path.exists(MEGAHIT):
  os.makedirs(MEGAHIT)

BENCH = os.path.join(OUTPUT, 'BENCHMARKS')
SUMDIR = os.path.join('hecatomb_report')
ASSEMBLY = os.path.join(WORKDIR, 'ASSEMBLY')
MAPPING = os.path.join(WORKDIR, 'MAPPING')
STATS = os.path.join(WORKDIR, 'STATS')
