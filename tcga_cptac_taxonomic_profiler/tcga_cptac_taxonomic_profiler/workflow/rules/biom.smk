##### add biom for the bracken reports next - script with glob



# kraken-biom *kraken2_bracken -o bracken_species.biom --fmt json


rule biom:
    input:
        expand(os.path.join(BRACKEN,"{sample}.kraken_bracken_genus.txt"), sample = SAMPLES),
        expand(os.path.join(BRACKEN,"{sample}.kraken_bracken_species.txt"), sample = SAMPLES)
    output:
        os.path.join(BIOM,"bracken_genus.biom"),
        os.path.join(BIOM,"bracken_species.biom")
    conda:
        os.path.join('..', 'envs','biom.yaml')
    log: 
        os.path.join(LOGS, "biom", "biom.log")
    benchmark: 
        os.path.join(BENCHMARKS, "biom", "biom.log")
    resources:
        mem_mb=SmallJobMem,
        time=60
    threads:
        1
    shell:
        '''
        kraken-biom {input[0]} -o {output[0]}  --fmt json
        kraken-biom {input[1]} -o {output[1]}  --fmt json
        '''

rule aggr_biom:
    """Aggregate biom"""
    input:
        os.path.join(BIOM,"bracken_genus.biom"),
        os.path.join(BIOM,"bracken_species.biom")
    output:
        os.path.join(FLAGS, "aggr_biom.flag")
    resources:
        mem_mb=SmallJobMem,
        time=5
    threads:
        1
    shell:
        """
        touch {output[0]}
        """