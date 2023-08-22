rule bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        os.path.join(BAMS_DIR, "{sample}.bam")
    output:
        os.path.join(BAMS_DIR,"{sample}.bam.bai")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        samtools index -@ {threads} {input[0]} {output[0]} 
        """

rule bam_unmap_sort_fastq:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(BAMS_DIR,"{sample}.bam")
    output:
        os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"),
        os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    shell:
        """
        samtools view -u -f 12 -F 256 -@ {threads} {input[0]} | samtools sort -@ {threads} |   
        samtools fastq -@ {threads} \
        -1 {output[0]} \
        -2 {output[1]} \
        -0 /dev/null -s /dev/null -n 
        """



#### aggregation rule

rule aggr_bam_to_fastq:
    """Aggregate for flag."""
    input:
        expand(os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_fastq.flag")
    threads:
        1
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """
