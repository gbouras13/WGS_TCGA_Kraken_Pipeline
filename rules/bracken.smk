rule bracken:
    input:
        os.path.join(RESULTS,"{sample}.kraken.txt"),
        os.path.join(RESULTS,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    output:
        os.path.join(RESULTS,"{sample}.kraken_bracken_genuses.txt"),
        os.path.join(RESULTS,"{sample}.kraken_bracken_species.txt")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    shell:
        '''
        bracken -d {params[0]} -i {input[1]} -o {output[1]}  -r 50 -l S 
        bracken -d {params[0]} -i {input[1]}  -o {output[0]} -r 50 -l G 
        '''
# https://github.com/jenniferlu717/Bracken/issues/81
# use 50


rule biom:
    input:
        expand(os.path.join(RESULTS,"{sample}.kraken_bracken_genuses.txt"), sample = SAMPLES),
        expand(os.path.join(RESULTS,"{sample}.kraken_bracken_species.txt"), sample = SAMPLES)
    output:
        os.path.join(BIOM,"bracken_genus.biom"),
        os.path.join(BIOM,"bracken_species.biom")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    shell:
        '''
        kraken-biom {input[0]} -o {output[0]}  --fmt json
        kraken-biom {input[1]}  -o {output[1]} --fmt json
        '''


rule aggr_bracken:
    """aggregated"""
    input:
        os.path.join(BIOM,"bracken_species.biom")
    output:
        os.path.join(LOGS, "aggr_bracken.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """