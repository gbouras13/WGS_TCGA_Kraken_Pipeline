rule bracken_s:
    input:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    output:
        os.path.join(KRAKEN_S,"{sample}.kraken_bracken_species.txt")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    shell:
        '''
        bracken -d {params[0]} -i {input[1]} -o {output[0]}  -r 50 -l S 
        '''
# https://github.com/jenniferlu717/Bracken/issues/81
# use 50

# rule bracken_g:
#     input:
#         os.path.join(KRAKEN_G,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_G,"{sample}.kraken.rep")
#     params:
#         os.path.join(DBDIR, 'standard')
#     output:
#         os.path.join(KRAKEN_G,"{sample}.kraken_bracken_genuses.txt")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     shell:
#         '''
#         bracken -d {params[0]} -i {input[1]}  -o {output[0]} -r 50 -l G 
#         '''

##### add biom for the bracken reports next - script with glob



# combine_bracken_outputs.py --files kraken/*total.txt -o bracken_species_all
# kraken-biom bracken/*_report_bracken_species.txt -o bracken_species.biom --fmt json

##### do this 

# rule combine_bracken_s:
#     """combine"""
#     input:
#         expand(os.path.join(KRAKEN_S,"{sample}.kraken_bracken_species.txt"), sample = SAMPLES)
#     output:
#         os.path.join(BRACKEN, "aggr_bracken.txt")
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         touch {output[0]}
#         """

# combine_bracken_outputs.py --files kraken/*total.txt -o bracken_species_all
# kraken-biom bracken/*_report_bracken_species.txt -o bracken_species.biom --fmt json

# rule biom:
#     input:
#         expand(os.path.join(RESULTS,"{sample}.kraken_bracken_genuses.txt"), sample = SAMPLES),
#         expand(os.path.join(RESULTS,"{sample}.kraken_bracken_species.txt"), sample = SAMPLES)
#     output:
#         os.path.join(BIOM,"bracken_genus.biom"),
#         os.path.join(BIOM,"bracken_species.biom")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     shell:
#         '''
#         kraken-biom {input[0]} -o {output[0]}  --fmt json
#         kraken-biom {input[1]}  -o {output[1]} --fmt json
#         '''


rule aggr_bracken:
    """aggregated"""
    input:
        expand(os.path.join(KRAKEN_S,"{sample}.kraken_bracken_species.txt"), sample = SAMPLES)
        # ,
        # expand(os.path.join(KRAKEN_S,"{sample}.kraken_bracken_species.txt"), sample = SAMPLES)
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