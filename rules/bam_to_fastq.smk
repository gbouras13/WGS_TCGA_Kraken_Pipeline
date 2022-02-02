
#### if you want to only look at the unmapped READS
#### I have decided to look at them all for now

rule bam_index:
    """Index a .bam file for rapid access with samtools."""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(READS,"{sample}.bam.bai")
    log:
        os.path.join(RESULTS,"{sample}.samtools.stderr")
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
        os.path.join(RESULTS,"{sample}.bam_unmap.log")
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

rule bam_sort:
    """Sort Unmapped Reads"""
    input:
        os.path.join(TMP,"{sample}_unmap.bam")
    output:
        os.path.join(TMP,"{sample}_unmap_sorted.bam")
    log:
        os.path.join(RESULTS,"{sample}.bam_sort.log")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools sort -@ {threads} {input[0]} > {output[0]} 2> {log}
        """

rule bam_to_fastq:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    log:
        os.path.join(RESULTS,"{sample}.bam_to_fastq.log")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools fastq -@ {threads} {input[0]} \
        -1 {output[0]} \
        -2 {output[1]} \
        -0 /dev/null -s /dev/null -n 2> {log}
        """


rule kraken2:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    log:
        os.path.join(RESULTS,"{sample}.bam_to_fastq.log")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools fastq -@ {threads} {input[0]} \
        -1 {output[0]} \
        -2 {output[1]} \
        -0 /dev/null -s /dev/null -n 2> {log}
        """

#### aggregation rule

rule test:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(TMP,"{sample}_R1.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(RESULTS, "test.txt")
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
