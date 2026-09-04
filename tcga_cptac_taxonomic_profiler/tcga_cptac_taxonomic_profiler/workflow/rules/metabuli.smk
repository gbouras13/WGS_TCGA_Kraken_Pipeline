rule run_metabuli:
    """Classify host-depleted reads with Metabuli.

    Metabuli is algorithmically distinct from Kraken2: it indexes 'metamers'
    spanning both amino acid and DNA space, using AA conservation for sensitive
    homology detection and DNA substitutions to separate close relatives.
    Kraken2 matches exact DNA k-mers only.

    This matters most for the 51 bp samples, which are 60% of the primary
    tumours and where exact DNA k-mer matching is weakest -- precisely the
    regime that generates the false positives the cancer-microbiome literature
    was retracted over.

    --precise 1 selects the short-read preset. The report is emitted in Kraken2
    format, so it feeds the existing report-parsing path unchanged.

    CLI (confirmed against metabuli 1.2.0):
      metabuli classify <query file(s)> <DB dir> <out dir> <job ID> [options]
    --seq-mode 2 is paired-end. It is the default, but set explicitly so the
    run does not depend on an upstream default staying put.
    """
    input:
        r1 = os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz"),
        r2 = os.path.join(INPUT, "{sample}_R2.host_rm.fastq.gz")
    output:
        report = os.path.join(METABULI, "{sample}_report.tsv"),
        classif = os.path.join(METABULI, "{sample}_classifications.tsv")
    params:
        db = config.databases.metabuli,
        outdir = METABULI,
        job_id = "{sample}",
        max_ram = config.metabuli.max_ram_gib,
        precise = config.metabuli.precise,
        min_score = config.metabuli.min_score,
        seq_mode = config.metabuli.seq_mode
    conda:
        os.path.join('..', 'envs', 'metabuli.yaml')
    log:
        os.path.join(LOGS, "metabuli", "{sample}.metabuli.log")
    benchmark:
        os.path.join(BENCHMARKS, "metabuli", "{sample}.metabuli.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    shell:
        """
        metabuli classify \
            {input.r1} {input.r2} \
            {params.db} \
            {params.outdir} \
            {params.job_id} \
            --threads {threads} \
            --seq-mode {params.seq_mode} \
            --max-ram {params.max_ram} \
            --precise {params.precise} \
            --min-score {params.min_score} 2> {log}
        """


rule aggr_metabuli:
    """Aggregate metabuli."""
    input:
        expand(os.path.join(METABULI, "{sample}_report.tsv"), sample = SAMPLES),
        expand(os.path.join(METABULI, "{sample}_classifications.tsv"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_metabuli.flag")
    threads:
        1
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """
