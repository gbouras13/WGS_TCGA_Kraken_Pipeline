rule run_kraken_s:
    """Index a .bam file for rapid access with samtools."""
    input:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    log:
        os.path.join(LOGS,"{sample}.kraken_s.log")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        kraken2 {input[0]} {input[1]}  \
                --threads {threads} --db {params[0]} --output {output[0]} \
                --paired \
                --report-minimizer-data \
                --confidence 0.05 --report {output[1]} 2> {log}
        """

# rule run_kraken_g:
#     """Index a .bam file for rapid access with samtools."""
#     input:
#         os.path.join(TMP,"{sample}_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_R2.fastq.gz")
#     output:
#         os.path.join(KRAKEN_G,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_G,"{sample}.kraken.rep")
#     params:
#         os.path.join(DBDIR, 'standard')
#     log:
#         os.path.join(LOGS,"{sample}.kraken_g.log")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         kraken2 {input[0]} {input[1]}  \
#                 --threads {threads} --db {params[0]} --output {output[0]} \
#                 --paired \
#                 --report-minimizer-data \
#                 --confidence 0.05 --report {output[1]} 2> {log}
#         """

# python3 ../../WGS_TCGA_Kraken_Pipeline/Kraken_Tools/extract_kraken_reads.py -k TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted.kraken.txt -s1 ../PROCESSING/TMP/TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted_R1.fastq.gz -s2 ../PROCESSING/TMP/TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted_R2.fastq.gz -o test_out.fastq -r TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted.kraken.rep  -t 2 --include-children  --fastq-output

# python3 ../../WGS_TCGA_Kraken_Pipeline/Kraken_Tools/extract_kraken_reads.py -k TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted.kraken.txt -s1 ../PROCESSING/TMP/TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted_R1.fastq.gz -s2 ../PROCESSING/TMP/TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted_R2.fastq.gz -o test_out_R1.fastq -o2 test_out_R2.fastq -r TCGA-CN-5361-01A-01D-1431_120418_SN1222_0098_BD0T68ACXX_s_2_rg.sorted.kraken.rep  -t 2 --include-children  --fastq-output

# python3 ../../WGS_TCGA_Kraken_Pipeline/Kraken_Tools/extract_kraken_reads.py -k TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted.kraken.txt -s1 TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted_R1.fastq.gz -s2 TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted_R2.fastq.gz -o test_out_R1.fastq -o2 test_out_R2.fastq -r TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted.kraken.rep  -t 2 --include-children  --fastq-output


# megahit -1 test_out_R1.fastq -2 test_out_R2.fastq -o out

# gondiiii

# python3 ../../WGS_TCGA_Kraken_Pipeline/Kraken_Tools/extract_kraken_reads.py -k TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted.kraken.txt -s1 TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted_R1.fastq.gz -s2 TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted_R2.fastq.gz -o test_out_R1_gondi.fastq -o2 test_out_R2_gondi.fastq -r TCGA-BA-5149-01A-01D-1509_120420_SN1120_0135_BD0T5FACXX_s_1_rg.sorted.kraken.rep -t 5811 --include-children  --fastq-output

# megahit -1 test_out_R1_gondi.fastq  -2 test_out_R2_gondi.fastq  -o gondi



rule aggr_kraken:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(KRAKEN_S,"{sample}.kraken.txt"), sample = SAMPLES),
        expand(os.path.join(KRAKEN_S,"{sample}.kraken.rep"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_kraken.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
