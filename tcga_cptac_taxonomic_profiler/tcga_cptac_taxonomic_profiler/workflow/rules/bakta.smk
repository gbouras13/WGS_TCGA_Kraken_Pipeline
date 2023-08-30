rule run_bakta:
    """
    generate mash db
    """
    input:
        mag =  os.path.join(ALL_MAGS, '{mag}.fna')
    output:
        out_tsv = os.path.join(BAKTA, '{mag}', '{mag}.tsv') 
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    conda:
        os.path.join("..", "envs", "bakta.yaml")
    benchmark:
        os.path.join(BENCHMARKS, 'bakta', "{mag}.txt")
    log:
        os.path.join(LOGS, 'bakta', "{mag}.log")
    params:
        outdir = os.path.join(BAKTA, '{mag}'),
        db=config.databases.bakta,
    shell:
        """
        bakta --db {params.db} --output {params.outdir} -f -t {threads}  {input.mag}
        """


rule aggr_bakta:
    """
    generate mash db
    """
    input:
        tsvs = expand(os.path.join(BAKTA, '{mag}', '{mag}.tsv') , mag=HQ_MED_MAGS)
    output:
        outtouch = os.path.join(FLAGS, 'bakta.flag'),
    threads:
        config.resources.sml.cpu
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output.outtouch}

        """
