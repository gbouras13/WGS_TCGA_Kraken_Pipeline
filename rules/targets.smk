"""
All target output files for Hecatomb are declared here
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
    #os.path.join(LOGS, "aggr_kraken_second_pass.txt"),
    #os.path.join(LOGS, "aggr_bracken.txt"),

]

# if skipAssembly:
#     AssemblyFiles = []
#     ContigAnnotFiles = []
#     MappingFiles = []
# else:
#     # Assembly files
#     AssemblyFiles = [
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","assembly.fasta"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING","contig_count_table.tsv"),
#         os.path.join(RESULTS,"assembly.properties.tsv")]
#     # Contig annotations
#     ContigAnnotFiles = [
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_phylum_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_class_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_order_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_family_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_genus_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_species_summary.tsv")
#     ]
#     # Mapping files
#     MappingFiles = [
#         os.path.join(MAPPING,"assembly.seqtable.bam"),
#         os.path.join(MAPPING,"assembly.seqtable.bam.bai"),
#         os.path.join(RESULTS,"contigSeqTable.tsv"),
#         os.path.join(SUMDIR,"contigKrona.html")
#     ]
#
#
# # Secondary AA search files
# SecondarySearchFilesAA = [
#     os.path.join(SECONDARY_AA_OUT, "AA_bigtable.tsv"),
# ]
#
# # Secondary NT search files
# SecondarySearchFilesNT = [
#     os.path.join(SECONDARY_NT_OUT, "NT_bigtable.tsv"),
#     os.path.join(RESULTS, "bigtable.tsv"),
# ]
#
# # Summary files
# optionalSummary = [
#     os.path.join(SUMDIR,"Step00_counts.tsv"),
#     os.path.join(SUMDIR,"Step01_counts.tsv"),
#     os.path.join(SUMDIR,"Step02_counts.tsv"),
#     os.path.join(SUMDIR,"Step03_counts.tsv"),
#     os.path.join(SUMDIR,"Step04_counts.tsv"),
#     os.path.join(SUMDIR,"Step05_counts.tsv"),
#     os.path.join(SUMDIR,"Step06_counts.tsv"),
#     os.path.join(SUMDIR,"Step07_counts.tsv"),
#     os.path.join(SUMDIR,"Step08_counts.tsv"),
#     os.path.join(SUMDIR,"Step09_counts.tsv"),
#     os.path.join(SUMDIR,"Step10_counts.tsv"),
#     os.path.join(SUMDIR,"Step11_counts.tsv"),
#     os.path.join(SUMDIR,"Step12_counts.tsv"),
#     os.path.join(SUMDIR,"Step13_counts.tsv"),
#     os.path.join(SUMDIR,"Sankey.svg"),
# ]
# SummaryFiles = [
#     optionalSummary,
#     os.path.join(SUMDIR, 'hecatomb.samples.tsv'),
#     os.path.join(SUMDIR, "taxonLevelCounts.tsv"),
#     os.path.join(SUMDIR, "krona.html")
# ]
