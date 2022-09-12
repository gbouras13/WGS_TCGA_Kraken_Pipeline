"""
All target output files are declared here
"""

# Preprocessing files
PreprocessingFiles = [
    # expand(os.path.join(KRAKEN_S,"{sample}.kraken.txt"), sample = SAMPLES),
    # expand(os.path.join(KRAKEN_S,"{sample}.kraken.rep"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_R1.fastq.gz"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_R2.fastq.gz"), sample = SAMPLES),
    # expand(os.path.join(KRAKEN_S,"{sample}.kraken_second_pass.txt"), sample = SAMPLES),
    # expand(os.path.join(KRAKEN_S,"{sample}.kraken_second_pass.rep"), sample = SAMPLES),
    # expand(os.path.join(KRAKEN_S,"{sample}.kraken_bracken_species_second_pass.txt"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_bacteria_R1.fastq"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_bacteria_R2.fastq"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_virus_R1.fastq"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_virus_R2.fastq"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_bacteria_fastp_R1.fastq.gz"), sample = SAMPLES),
    # expand(os.path.join(TMP,"{sample}_bacteria_fastp_R2.fastq.gz"), sample = SAMPLES)


    os.path.join(LOGS, "aggr_fastq.txt"),
    os.path.join(LOGS, "aggr_kraken.txt"),
    os.path.join(LOGS, "aggr_assembly.txt"),
    os.path.join(LOGS, "aggr_fastp.txt"),
    os.path.join(LOGS, "aggr_kraken_second_pass.txt"),
    os.path.join(LOGS, "aggr_bracken.txt"),
    os.path.join(LOGS, "aggr_read_count.txt")

]
