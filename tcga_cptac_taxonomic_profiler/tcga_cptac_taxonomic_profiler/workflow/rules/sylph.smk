rule run_sylph:
    """Profile host-depleted reads with sylph (ANI-based containment).

    Sylph is the fourth profiler and answers a question the others cannot:
    *how similar is the detected genome to its reference*. Kraken2 and Metabuli
    report that a taxon is present; sylph reports the containment ANI, so a hit
    at 99% ANI is the named organism while a hit at 90% is k-mer noise wearing a
    species name. That distinction is the substance of the cancer-microbiome
    retractions.

    It is also built for exactly this regime: zero-inflated Poisson debiasing of
    ANI at low coverage. Median load here is 0.06-0.33%, and SingleM returned
    nothing for 196 of 336 samples.

    CLI verified against sylph 0.9.0:
      sylph profile <db.syldb> -1 R1 -2 R2 -t N -o out.tsv

    --minimum-ani is set BELOW the 95 default on purpose. The low-ANI tail is
    the contamination signal we want to observe and filter in R, not something
    to discard silently at profiling time.
    """
    input:
        r1 = os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz"),
        r2 = os.path.join(INPUT, "{sample}_R2.host_rm.fastq.gz"),
        db = config.databases.sylph
    output:
        tsv = os.path.join(SYLPH, "{sample}.sylph.tsv")
    params:
        min_ani = config.sylph.min_ani,
        min_kmers = config.sylph.min_kmers,
        read_counts = "--estimate-read-counts" if config.sylph.estimate_read_counts else "",
        unknown = "-u" if config.sylph.estimate_unknown else ""
    conda:
        os.path.join('..', 'envs', 'sylph.yaml')
    log:
        os.path.join(LOGS, "sylph", "{sample}.sylph.log")
    benchmark:
        os.path.join(BENCHMARKS, "sylph", "{sample}.sylph.log")
    resources:
        # sylph is far lighter than kraken2 (~30x less memory in the paper)
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.med.cpu
    shell:
        """
        sylph profile {input.db} \
            -1 {input.r1} -2 {input.r2} \
            -t {threads} \
            --minimum-ani {params.min_ani} \
            --min-number-kmers {params.min_kmers} \
            {params.read_counts} {params.unknown} \
            -o {output.tsv} 2> {log}

        # sylph writes only a header when nothing passes threshold; keep the
        # file so the DAG completes, as the singlem rule does.
        touch {output.tsv}
        """


rule aggr_sylph:
    """Aggregate sylph."""
    input:
        expand(os.path.join(SYLPH, "{sample}.sylph.tsv"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_sylph.flag")
    threads:
        1
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """
