rule bam_counts_all:
    """Get read count from bam file."""
    input:
        os.path.join(BAMS_DIR, "{sample}.bam")
    output:
        os.path.join(READCOUNT,"{sample}_readcount_all.txt")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpus
    shell:
        """
        samtools view -c -@ {threads} {input[0]} > {output[0]}
        """


rule aggr_read_counts:
    """Aggregates Read Count."""
    input:
        expand(os.path.join(READCOUNT,"{sample}_readcount_all.txt"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_read_count.flag")
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
