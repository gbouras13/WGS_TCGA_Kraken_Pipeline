#### if you want to only look at the unmapped READS

rule bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        os.path.join(BAMS, "{sample}.bam")
    output:
        os.path.join(BAMS,"{sample}.bam.bai")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem,
        time=60
    shell:
        """
        samtools index -@ {threads} {input[0]} {output[0]} 
        """

rule bam_unmap_sort_fastq:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(BAMS,"{sample}.bam")
    output:
        os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"),
        os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
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
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"), sample = SAMPLES)
        expand(os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_fastq.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=5
    shell:
        """
        touch {output[0]}
        """
