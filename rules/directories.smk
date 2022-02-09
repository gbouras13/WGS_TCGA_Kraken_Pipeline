"""
Database and output locations for Hecatomb
Ensures consistent variable names and file locations for the pipeline, the database download script,
and the addHost script.
"""


DBDIR = 'Databases'

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


### MMSEQS DBs
UNIVIRDB = os.path.join(DBDIR, "aa", "virus_primary_aa")
UNIREF50VIR = os.path.join(DBDIR, "aa", "virus_secondary_aa")
NCBIVIRDB = os.path.join(DBDIR, "nt", "virus_primary_nt")
POLYMICRODB = os.path.join(DBDIR, "nt", "virus_secondary_nt")


### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')
WORKDIR = os.path.join(OUTPUT, 'PROCESSING')
TMP = os.path.join(WORKDIR, 'TMP')
LOGS = os.path.join(OUTPUT, 'LOGS')
BIOM = LOGS = os.path.join(RESULTS, 'BIOM')



BENCH = os.path.join(OUTPUT, 'BENCHMARKS')
SUMDIR = os.path.join('hecatomb_report')
ASSEMBLY = os.path.join(WORKDIR, 'ASSEMBLY')
MAPPING = os.path.join(WORKDIR, 'MAPPING')
STATS = os.path.join(WORKDIR, 'STATS')
