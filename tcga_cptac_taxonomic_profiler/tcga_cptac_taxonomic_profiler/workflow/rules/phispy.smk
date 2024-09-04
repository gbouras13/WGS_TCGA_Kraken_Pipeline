
rule phispy:
    """Run phispy."""
    input:
        gbk = os.path.join(BAKTA, '{mag}', '{mag}.gbff') 
    output:
        outcoord = os.path.join(PHISPY,"{mag}", "prophage_coordinates.tsv"),
        outfasta = os.path.join(PHISPY,"{mag}", "phage.fasta")
    params:
        outdir = os.path.join(PHISPY,"{mag}")
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    conda:
        os.path.join("..", "envs", "phispy.yaml")
    benchmark:
        os.path.join(BENCHMARKS, 'phispy', "{mag}.txt")
    log:
        os.path.join(LOGS, 'phispy', "{mag}.log")
    shell:
        """
        phispy {input.gbk} --output_choice 512 -o {params.outdir} --phage_genes 0
        """


rule aggr_phispy:
    """Aggregate."""
    input:
        expand(os.path.join(PHISPY,"{mag}", "prophage_coordinates.tsv"), mag = HQ_MED_MAGS),
        expand(os.path.join(PHISPY,"{mag}", "phage.fasta"), mag = HQ_MED_MAGS)
    output:
        outtouch = os.path.join(FLAGS, 'phispy.flag'),
    threads:
        config.resources.sml.cpu
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """
