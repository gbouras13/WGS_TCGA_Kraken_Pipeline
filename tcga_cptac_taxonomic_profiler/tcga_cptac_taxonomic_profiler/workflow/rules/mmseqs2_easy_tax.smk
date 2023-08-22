
"""

Taken from Rob's [atavide_lite](https://github.com/linsalrob/atavide_lite/blob/main/slurm/mmseqs_easy_taxonomy.slurm)

"""

rule run_mmseqs_easy_tax:
    input: 
        fasta1 = os.path.join(FASTA, "{sample}.R1.fasta"),
        fasta2 = os.path.join(FASTA, "{sample}.R2.fasta"),
    output:
        outtouch=os.path.join(MMSEQS2, 'flags', '{sample}.done'),
        tmp = temp(os.path.join(dir.out.assembly, "coAssembly", "co_assembly_graph.fastg")),
        tmpdir = temp(os.path.join(TMPDIR, "{sample}"))
    params:
        db = config.databases.mmseqs2.uniref50,
        outdir=os.path.join(MMSEQS2, '{sample}')
    threads: 
        config.resources.big.cpu
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.big.time
    log: 
        os.path.join(LOGS, "mmseqs2", "{sample}.taxonomy.log")
    benchmark: 
        os.path.join(BENCHMARKS, "mmseqs2", "{sample}.taxonomy.benchmark")
    conda: 
        os.path.join("..", "envs", "mmseqs2.yaml")
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        mmseqs easy-taxonomy {input.fasta1} {input.fasta2} {params.db} {params.outdir} {output.tmpdir} --start-sens 1 --sens-steps 3 -s 7 --threads {threads} --orf-filter 0 2>> {log}
        touch {output.outtouch}
        """


rule aggr_mmseqs2_easy_tax:
    input:
        expand(os.path.join(MMSEQS2, 'flags', '{sample}.done'), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_mmseqs2.flag")
    threads:
        1
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """