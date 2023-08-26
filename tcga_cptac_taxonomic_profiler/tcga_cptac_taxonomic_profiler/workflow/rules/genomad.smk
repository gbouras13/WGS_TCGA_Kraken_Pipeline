"""
use genomad
must be run after the sample_assembly_binning step (as it runs it on the whole catalogue of contigs)
"""

rule run_genomad:
    """
    https://github.com/apcamargo/genomad
    """

    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz')
    params:
        db = config.mge.genomad_db,
        outdir = GENOMAD_RESULTS
    output:
        outtouch = os.path.join(FLAGS, 'genomad.flag')
    benchmark:
        os.path.join(BENCHMARKS, 'genomad', "genomad.txt")
    log:
        os.path.join(LOGS, 'genomad', "genomad.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join("..", "envs", "genomad.yaml")
    shell:
        """
        genomad end-to-end --cleanup  {input.catalogue} {params.outdir} {params.db}
        touch {output.outtouch}
        """
