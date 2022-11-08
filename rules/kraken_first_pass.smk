
rule run_kraken_first_pass:
    """Runs kraken down to the species level.
    https://www.biorxiv.org/content/10.1101/2022.04.27.489753v1.full.pdf
    """
    input:
        os.path.join(UNALIGNED_FASTQ,"{sample}.R1.trimmed.fastq.gz"),
        os.path.join(UNALIGNED_FASTQ,"{sample}.R2.trimmed.fastq.gz")
    output:
        os.path.join(KRAKEN_FIRST_PASS,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_FIRST_PASS,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem,
        time=120
    shell:
        """
        kraken2 {input[0]} {input[1]}  \
                --threads {threads} --db {params[0]} --output {output[0]} \
                --paired \
                --report-minimizer-data \
                --confidence 0.15 --report {output[1]} 
        """
        

rule aggr_kraken_first_pass:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(KRAKEN_FIRST_PASS,"{sample}.kraken.txt"), sample = SAMPLES),
        expand(os.path.join(KRAKEN_FIRST_PASS,"{sample}.kraken.rep"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_kraken.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=5
    shell:
        """
        touch {output[0]}
        """
