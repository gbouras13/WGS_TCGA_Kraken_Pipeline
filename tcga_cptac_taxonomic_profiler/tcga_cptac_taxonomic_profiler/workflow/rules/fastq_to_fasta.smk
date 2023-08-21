rule fastq_to_fasta_R1:
    input:
        os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz")
    output:
        os.path.join(FASTA, "{sample}.R1.fasta")
    threads:
        1
    log: 
        os.path.join(LOGS, "fastq_to_fastq", "{sample}.log")
    benchmark: 
        os.path.join(BENCHMARKS, "fastq_to_fastq", "{sample}.benchmark")
    conda: 
        os.path.join("..", "envs", "biopython.yml")
    resources:
        mem_mb=SmallJobMem,
        time=MediumJobTimeMin
    script:
        '../scripts/fastq_to_fasta.py'

rule fastq_to_fasta_R2:
    input:
        os.path.join(INPUT, "{sample}_R2.host_rm.fastq.gz")
    output:
        os.path.join(FASTA, "{sample}.R2.fasta")
    threads:
        1
    log: 
        os.path.join(LOGS, "fastq_to_fastq", "{sample}.log")
    benchmark: 
        os.path.join(BENCHMARKS, "fastq_to_fastq", "{sample}.benchmark")
    conda: 
        os.path.join("..", "envs", "biopython.yml")
    resources:
        mem_mb=SmallJobMem,
        time=MediumJobTimeMin
    script:
        '../scripts/fastq_to_fasta.py'

rule aggr_fastq_to_fasta:
    input:
        expand(os.path.join(FASTA, "{sample}.fasta"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_fastq_to_fasta.flag")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=SmallJobTimeMin
    shell:
        """
        touch {output[0]}
        """
