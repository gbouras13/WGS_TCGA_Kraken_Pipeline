"""
All target output files are declared here
"""

# Preprocessing files
OutputFiles = [
    os.path.join(LOGS, "aggr_kraken.txt"),
    os.path.join(LOGS, "aggr_extraction.txt"),
    os.path.join(LOGS, "aggr_fastp.txt"),
    os.path.join(LOGS, "aggr_kraken_second_pass.txt"),
    os.path.join(LOGS, "aggr_bracken.txt")

]

# Preprocessing files
BamToFastqFiles = [
    os.path.join(LOGS, "aggr_fastq.txt"),
    os.path.join(LOGS, "aggr_read_count.txt")
]
