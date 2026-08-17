rule deacon_host_removal:
    """Deplete human reads with Deacon against the panhuman-1 pangenome index.

    panhuman-1 = HPRC year-1 assemblies + CHM13v2.0 + GRCh38.p14, with bacterial
    and viral sequence removed so genuine microbial reads are not depleted along
    with host. This supersedes depletion against a single linear reference: it
    additionally captures human haplotypes absent from any one assembly, which
    are a known source of false microbial classification.

    Produces the *.host_rm.fastq.gz files that the kraken, singlem, assembly and
    binning stages consume. Previously this step was performed outside the
    workflow, so its parameters were not recorded.

    Deacon defaults are -a 2 -r 0.01. We default to abs_threshold 1 because in
    tumour WGS the cost of retaining a human read (a false microbial call)
    greatly exceeds the cost of discarding a true microbial read. Run both and
    report the difference.
    """
    input:
        r1 = os.path.join(INPUT, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(INPUT, "{sample}_R2.fastq.gz"),
        index = os.path.join(config.databases.host_db, "panhuman-1.k31w15.idx")
    output:
        r1 = os.path.join(HOST_RM_FASTQ, "{sample}_R1.host_rm.fastq.gz"),
        r2 = os.path.join(HOST_RM_FASTQ, "{sample}_R2.host_rm.fastq.gz"),
        rs = os.path.join(HOST_RM_FASTQ, "{sample}_RS.host_rm.fastq.gz"),
        summary = os.path.join(HOST_RM_FASTQ, "{sample}.deacon.json")
    params:
        abs_threshold = config.deacon.abs_threshold,
        rel_threshold = config.deacon.rel_threshold
    conda:
        os.path.join('..', 'envs', 'deacon.yaml')
    log:
        os.path.join(LOGS, "deacon", "{sample}.deacon.log")
    benchmark:
        os.path.join(BENCHMARKS, "deacon", "{sample}.deacon.log")
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        deacon filter \
            --deplete \
            --threads {threads} \
            --abs-threshold {params.abs_threshold} \
            --rel-threshold {params.rel_threshold} \
            --summary {output.summary} \
            {input.index} {input.r1} {input.r2} \
            --output {output.r1} \
            --output2 {output.r2} 2> {log}

        # Deacon retains or discards read pairs together, so no singletons are
        # produced by construction. The extract stage also discards singletons
        # (samtools fastq -0 /dev/null -s /dev/null). An empty file is emitted
        # to satisfy the assembly stage, which expects an _RS input.
        : | gzip -c > {output.rs}
        """


rule aggr_host_removal:
    """Aggregate host removal."""
    input:
        expand(os.path.join(HOST_RM_FASTQ, "{sample}_R1.host_rm.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(HOST_RM_FASTQ, "{sample}_R2.host_rm.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(HOST_RM_FASTQ, "{sample}_RS.host_rm.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(HOST_RM_FASTQ, "{sample}.deacon.json"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_host_removal.flag")
    threads:
        1
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """
