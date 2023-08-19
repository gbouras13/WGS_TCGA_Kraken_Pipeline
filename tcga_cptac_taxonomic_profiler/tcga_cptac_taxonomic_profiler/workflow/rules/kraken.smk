
rule run_kraken_first_pass:
    """Runs kraken down to the species level.
    confidence based on this 
    https://www.biorxiv.org/content/10.1101/2022.04.27.489753v1.full.pdf
    """
    input:
        os.path.join(INPUT, "{sample}.R1.fastq.gz"),
        os.path.join(INPUT, "{sample}.R2.fastq.gz")
    output:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep")
    params:
        KRAKENDB
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    resources:
        mem_mb=BigJobMem,
        time=120
    threads:
        2
    shell:
        """
        kraken2 {input[0]} {input[1]}  \
                --threads {threads} --db {params[0]} --output {output[0]} \
                --paired \
                --report-minimizer-data \
                --confidence 0.15 --report {output[1]} 
        """
        

rule aggr_kraken:
    """Aggregate kraken"""
    input:
        expand(os.path.join(KRAKEN,"{sample}.kraken.txt"), sample = SAMPLES),
        expand(os.path.join(KRAKEN,"{sample}.kraken.rep"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_kraken.flag")
    resources:
        mem_mb=SmallJobMem,
        time=5
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
