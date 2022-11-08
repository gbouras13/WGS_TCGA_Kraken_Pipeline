rule run_kraken_s_second_pass:
    """run kraken on the fastq files that have been cleaned with fastp."""
    input:
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fastp_R1.fastq.gz"),
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fastp_R2.fastq.gz")
    output:
        os.path.join(KRAKEN_SECOND_PASS,"{sample}.kraken_second_pass.txt"),
        os.path.join(KRAKEN_SECOND_PASS,"{sample}.kraken_second_pass.rep")
    params:
        os.path.join(DBDIR, 'standard')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    resources:
        mem_mb=BigJobMem,
        time=120,
        th=BigJobCpu
    shell:
        """
        kraken2 {input[0]} {input[1]}  \
                --threads {resources.th} --db {params[0]} --output {output[0]} \
                --paired \
                --report-minimizer-data \
                --confidence 0.15 --report {output[1]} 
        """

rule aggr_kraken_second_pass:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(KRAKEN_SECOND_PASS,"{sample}.kraken_second_pass.txt"), sample = SAMPLES),
        expand(os.path.join(KRAKEN_SECOND_PASS,"{sample}.kraken_second_pass.rep"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_kraken_second_pass.txt")
    resources:
        mem_mb=SmallJobMem,
        time=5,
        th=1
    shell:
        """
        touch {output[0]}
        """
