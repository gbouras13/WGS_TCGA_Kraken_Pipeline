
rule run_singlem:
    """
    Runs singlem microbial fraction - note use the untrimmed FASTQs
    Need to touch the output files in case the job fails (aka not enough reads)
    """
    input:
        os.path.join(INPUT, "{sample}_R1.fastq.gz"),
        os.path.join(INPUT, "{sample}_R2.fastq.gz")
    output:
        os.path.join(SINGLEM,"{sample}.pipe.tsv"),
        os.path.join(SINGLEM,"{sample}.pipe.krona.txt"),
        os.path.join(SINGLEM,"{sample}.smf.tsv")
    params:
        config.databases.singlem
    conda:
        os.path.join('..', 'envs','singlem.yaml')
    log: 
        os.path.join(LOGS, "singlem", "{sample}.singlem.log")
    benchmark: 
        os.path.join(BENCHMARKS, "singlem", "{sample}.singlem.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        export SINGLEM_METAPACKAGE_PATH={params[0]}

        singlem pipe \
            -1 {input[0]} -2 {input[1]} -p {output[0]} --taxonomic-profile-krona {output[1]} --threads {threads}

        singlem microbial_fraction \
            --forward {input[0]} --reverse {input[1]}  -p {output[0]} > {output[2]}
        
        touch {output[0]}
        touch {output[1]}
        touch {output[2]}
        """

rule aggr_singlem:
    """Aggregate single"""
    input:
        expand(os.path.join(SINGLEM,"{sample}.smf.tsv"), sample = SAMPLES),
        expand(os.path.join(SINGLEM,"{sample}.smf.tsv"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "singlem.flag")
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output[0]}
        """
