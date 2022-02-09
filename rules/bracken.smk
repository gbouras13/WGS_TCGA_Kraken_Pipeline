rule bracken:
    input:
        os.path.join(RESULTS,"{sample}.kraken.txt"),
        os.path.join(RESULTS,"{sample}.kraken.rep")
    params:
        os.path.join(DBDIR, 'standard')
    output:
        os.path.join(RESULTS,"{sample}._report_bracken_species.txt"),
        os.path.join(RESULTS,"{sample}._report_bracken_genus.txt")
    conda:
        os.path.join('..', 'envs','bracken.yaml')
    shell:
        '''
        bracken -d {params.kraken_db} -i {input.report} -o {input.total} -r 50 -l S
        bracken -d {params.kraken_db} -i {input.report} -o {input.total} -r 50 -l G
        '''
# https://github.com/jenniferlu717/Bracken/issues/81
# use 50

rule aggr_bracken:
    """aggregated"""
    input:
        expand(os.path.join(RESULTS,"{sample}._report_bracken_species.txt"), sample = SAMPLES)
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
