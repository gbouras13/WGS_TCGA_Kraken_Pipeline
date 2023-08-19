# https://github.com/jenniferlu717/Bracken/issues/81
# use 50 as there is a mixture for consistency

rule bracken_species:
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep")
    params:
        KRAKENDB
    output:
        os.path.join(BRACKEN,"{sample}.kraken_bracken_species.txt")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    log: 
        os.path.join(LOGS, "bracken", "{sample}.bracken.species.log")
    benchmark: 
        os.path.join(BENCHMARKS, "bracken", "{sample}.bracken.species.log")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=60
    shell:
        '''
        bracken -d {params[0]} -i {input[1]} -o {output[0]}  -r 50 -l S 
        '''


rule bracken_genus:
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep")
    params:
        KRAKENDB
    output:
        os.path.join(BRACKEN,"{sample}.kraken_bracken_genus.txt")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    log: 
        os.path.join(LOGS, "bracken", "{sample}.bracken.genus.log")
    benchmark: 
        os.path.join(BENCHMARKS, "bracken", "{sample}.bracken.genus.log")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=60
    shell:
        '''
        bracken -d {params[0]} -i {input[1]}  -o {output[0]} -r 50 -l G 
        '''


rule aggr_bracken:
    """aggregated"""
    input:
        expand(os.path.join(BRACKEN,"{sample}.kraken_bracken_species.txt"), sample = SAMPLES),
        expand(os.path.join(BRACKEN,"{sample}.kraken_bracken_genus.txt"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_bracken.txt")
    resources:
        mem_mb=SmallJobMem,
        time=5
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
