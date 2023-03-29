rule host_removal_mapping:
    """Preprocessing step 02a: Host removal: mapping to host.
    
    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz"),
        host = HOSTINDEX
    output:
        r1 = os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R1.unmapped.fastq"),
        r2 = os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R2.unmapped.fastq"),
        s = os.path.join(HOST_UNMAPPED_FASTQ_SINGLETONS,  "{sample}_R1.unmapped.singletons.fastq"),
        o = os.path.join(HOST_UNMAPPED_FASTQ_SINGLETONS,  "{sample}_R1.other.singletons.fastq")
    resources:
        mem_mb = 16000,
        time = 300
    threads:
        16
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} --secondary=no {input.host} {input.r1} {input.r2} \
            | samtools view -f 4 -h  \
            | samtools fastq -NO -1 {output.r1} -2 {output.r2} -0 {output.o} -s {output.s} 
        """


#### aggregation rule

rule aggr_host_removal:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R1.unmapped.fastq"), sample = SAMPLES),
        expand(os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R2.unmapped.fastq"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_host_removal.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=1
    shell:
        """
        touch {output[0]}
        """
