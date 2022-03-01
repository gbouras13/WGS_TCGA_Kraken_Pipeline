#### if you want to only look at the unmapped READS
#### I have decided to look at them all for now

rule bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(READS,"{sample}.bam.bai")
    log:
        os.path.join(LOGS,"{sample}.samtools.stderr")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools index -@ {threads} {input[0]} {output[0]} 2> {log}
        """
rule bam_unmap:
    """Get unmapped reads"""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(TMP,"{sample}_unmap.bam")
    log:
        os.path.join(LOGS,"{sample}.bam_unmap.log")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools view -u -f 12 -F 256 -@ {threads} {input[0]} > {output[0]} 2> {log}
        """

rule bam_unmap_sort_fastq:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    log:
        os.path.join(LOGS,"{sample}.bam_unmap_sort_fastq.log")
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
        -0 /dev/null -s /dev/null -n 2> {log}
        """


#### aggregation rule

rule test:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(TMP,"{sample}_R1.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(TMP,"{sample}_R2.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "test.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
