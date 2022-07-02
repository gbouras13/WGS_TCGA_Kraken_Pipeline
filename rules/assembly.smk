rule extract_bact_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_bacteria_R1.fastq"),
        os.path.join(TMP,"{sample}_bacteria_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    log:
        os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 2 --include-children  --fastq-output
        """


# rule extract_plasmid_fastqs:
#     """Extract Fastas."""
#     input:
#         os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
#         os.path.join(TMP,"{sample}_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_R2.fastq.gz")
#     output:
#         os.path.join(TMP,"{sample}_other_R1.fastq"),
#         os.path.join(TMP,"{sample}_other_R2.fastq")
#     params:
#         os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
#     log:
#         os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
#         -o {output[0]} -o2 {output[1]} -r {input[1]} -t 28384 --include-children  --fastq-output
#         """

# rule extract_virus_fastqs:
#     """Extract Fastas."""
#     input:
#         os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
#         os.path.join(TMP,"{sample}_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_R2.fastq.gz")
#     output:
#         os.path.join(TMP,"{sample}_virus_R1.fastq"),
#         os.path.join(TMP,"{sample}_virus_R2.fastq")
#     params:
#         os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
#     log:
#         os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
#         -o {output[0]} -o2 {output[1]} -r {input[1]} -t 10239 --include-children  --fastq-output
#         """




# rule fastp_trim:
#     """remove adapters etc """
#     input:
#         os.path.join(TMP,"{sample}_bacteria_R1.fastq"),
#         os.path.join(TMP,"{sample}_bacteria_R2.fastq")
#     output:
#         os.path.join(TMP,"{sample}_bacteria_trim_R1.fastq"),
#         os.path.join(TMP,"{sample}_bacteria_trim_R2.fastq")
#     log:
#         os.path.join(LOGS,"{sample}_fastp.log")
#     conda:
#         os.path.join('..', 'envs','assembly.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -w {threads}
#         """


# rule megahit:
#     """run megahit."""
#     input:
#         os.path.join(TMP,"{sample}_bacteria_trim_R1.fastq"),
#         os.path.join(TMP,"{sample}_bacteria_trim_R2.fastq")
#     output:
#         os.path.join(MEGAHIT, "{sample}", "final.contigs.fa")
#     params:
#         directory(os.path.join(MEGAHIT, '{sample}'))
#     log:
#         os.path.join(LOGS,"{sample}.megahit.log")
#     conda:
#         os.path.join('..', 'envs','assembly.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
# 	    rm -rf {params[0]}
#         megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
#         """



rule aggr_assembly:
    """aggr"""
    input:
        expand(os.path.join(TMP,"{sample}_bacteria_R1.fastq"), sample = SAMPLES),
        #expand(os.path.join(TMP,"{sample}_virus_R1.fastq"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_assembly.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem
    shell:
        """
        touch {output[0]}
        """
