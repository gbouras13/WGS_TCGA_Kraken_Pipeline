rule run_kraken_s_second_pass:
    """run kraken on the bacterial fastp."""
    input:
        os.path.join(TMP,"{sample}_bacteria_fastp_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_bacteria_fastp_R2.fastq.gz")
    output:
        os.path.join(KRAKEN_S,"{sample}.kraken_second_pass.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken_second_pass.rep")
    params:
        os.path.join(DBDIR, 'standard')
    log:
        os.path.join(LOGS,"{sample}.kraken_s_second_pass.log")
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
                --confidence 0.15 --report {output[1]} 2> {log}
        """


rule aggr_kraken_second_pass:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(KRAKEN_S,"{sample}.kraken_second_pass.txt"), sample = SAMPLES),
        expand(os.path.join(KRAKEN_S,"{sample}.kraken_second_pass.rep"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_kraken_second_pass.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem
    shell:
        """
        touch {output[0]}
        """
