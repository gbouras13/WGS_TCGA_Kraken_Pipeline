rule bracken_species:
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    output:
        os.path.join(BRACKEN,"{sample}.kraken_bracken_species.txt")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=60
    shell:
        '''
        bracken -d {params[0]} -i {input[1]} -o {output[0]}  -r 50 -l S 
        '''
# https://github.com/jenniferlu717/Bracken/issues/81
# use 50

rule bracken_first_pass_genus:
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    output:
        os.path.join(BRACKEN,"{sample}.kraken_bracken_genus.txt")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=60
    shell:
        '''
        bracken -d {params[0]} -i {input[1]}  -o {output[0]} -r 50 -l G 
        '''

##### add biom for the bracken reports next - script with glob

# combine_bracken_outputs.py --files kraken/*total.txt -o bracken_species_all
# kraken-biom bracken/*_report_bracken_species.rep -o bracken_species.biom --fmt json

# kraken-biom *kraken2_bracken -o bracken_species.biom --fmt json

##### do this 
# combine_bracken_outputs.py --files kraken/*total.txt -o bracken_species_all
# kraken-biom bracken/*_report_bracken_species.txt -o bracken_species.biom --fmt json

# rule biom:
#     input:
#         expand(os.path.join(BRACKEN_FIRST_PASS,"{sample}.kraken_bracken_genus_first_pass.txt"), sample = SAMPLES),
#         expand(os.path.join(BRACKEN_FIRST_PASS,"{sample}.kraken_bracken_species_first_pass.txt"), sample = SAMPLES)
#     output:
#         os.path.join(BIOM,"bracken_genus_first_pass.biom"),
#         os.path.join(BIOM,"bracken_species_first_pass.biom")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     resources:
#         mem_mb=SmallJobMem,
#         time=60
#     threads:
#         1
#     shell:
#         '''
#         kraken-biom {input[0]} -o {output[0]}  --fmt json
#         kraken-biom {input[1]} -o {output[1]}  --fmt json
#         '''


rule aggr_bracken:
    """aggregated"""
    input:
        expand(os.path.join(BRACKEN,"{sample}.kraken_bracken_species_first_pass.txt"), sample = SAMPLES),
        expand(os.path.join(BRACKEN,"{sample}.kraken_bracken_genus_first_pass.txt"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_bracken.txt")
    resources:
        mem_mb=SmallJobMem,
        time=5
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
