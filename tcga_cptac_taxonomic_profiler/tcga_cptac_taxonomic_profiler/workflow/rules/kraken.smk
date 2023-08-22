
rule run_kraken:
    """Runs kraken down to the species level.
    confidence based on this 
    https://www.biorxiv.org/content/10.1101/2022.04.27.489753v1.full.pdf
    """
    input:
        os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz"),
        os.path.join(INPUT, "{sample}_R2.host_rm.fastq.gz")
    output:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep")
    params:
        config.databases.kraken
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    log: 
        os.path.join(LOGS, "kraken", "{sample}.kraken.log")
    benchmark: 
        os.path.join(BENCHMARKS, "kraken", "{sample}.kraken.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.med.time
    threads:
        config.resources.sml.cpu
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
        os.path.join(FLAGS, "aggr_kraken.flag")
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output[0]}
        """
