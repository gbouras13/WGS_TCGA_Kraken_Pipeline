rule bam_counts_all:
    """Get read count from bam file."""
    input:
        os.path.join(BAMS, "{sample}.bam")
    output:
        os.path.join(READCOUNT,"{sample}_readcount_all.txt")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    resources:
        mem_mb=MediumJobMem,
        time=600,
        th=16
    shell:
        """
        samtools view -c -@ {resources.th} {input[0]} > {output[0]}
        """


rule aggr_read_counts:
    """Aggregates Read Count."""
    input:
        expand(os.path.join(READCOUNT,"{sample}_readcount_all.txt"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_read_count.txt")
    resources:
        mem_mb=SmallJobMem,
        time=5,
        th=1
    shell:
        """
        touch {output[0]}
        """
