"""
Database and output locations for Hecatomb
Ensures consistent variable names and file locations for the pipeline, the database download script,
and the addHost script.
"""


DBDIR = 'Databases'
KRAKENTOOLSDIR = 'Kraken_Tools'
TRUST4DIR = 'Trust_4_Files'
CONTAMINANTS = 'contaminants'

### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'wgs_tcga_out'
else:
    OUTPUT = config['Output']


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
READCOUNT = = os.path.join(RESULTS, 'READCOUNT')



# needs to be created before megahit is run
if not os.path.exists(MEGAHIT):
  os.makedirs(MEGAHIT)

