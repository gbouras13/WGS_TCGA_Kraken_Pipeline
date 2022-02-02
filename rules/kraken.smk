rule run_kraken:
    """Index a .bam file for rapid access with samtools."""
    input:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(RESULTS,"{sample}.kraken.txt"),
        os.path.join(RESULTS,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    log:
        os.path.join(RESULTS,"{sample}.kraken.stderr")
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


rule test2:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(RESULTS,"{sample}.kraken.txt"), sample = SAMPLES)
    output:
        os.path.join(RESULTS, "test_2.txt")
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
