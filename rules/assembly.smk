rule extract_bact_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_bacteria_R1.fastq"),
        os.path.join(TMP,"{sample}_bacteria_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    log:
        os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 2 --include-children  --fastq-output
        """


rule extract_virus_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_virus_R1.fastq"),
        os.path.join(TMP,"{sample}_virus_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    log:
        os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 10239 --include-children  --fastq-output
        """


rule megahit:
    """run megahit."""
    input:
        os.path.join(TMP,"{sample}_bacteria_virus_fastp_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_bacteria_virus_fastp_R2.fastq.gz")
    output:
        os.path.join(MEGAHIT, "{sample}", "final.contigs.fa")
    params:
        directory(os.path.join(MEGAHIT, '{sample}'))
    log:
        os.path.join(LOGS,"{sample}.megahit.log")
    conda:
        os.path.join('..', 'envs','assembly.yaml')
    threads:
        8
    resources:
        mem_mb=MediumJobMem,
        time=60
    shell:
        """
	    rm -rf {params[0]}
        megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
        """


rule aggr_assembly:
    """aggr"""
    input:
        expand(os.path.join(TMP,"{sample}_bacteria_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(TMP,"{sample}_virus_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(MEGAHIT, "{sample}", "final.contigs.fa"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_assembly.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem
    shell:
        """
        touch {output[0]}
        """
