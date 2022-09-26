rule extract_bact_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
        os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"),
        os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz")
    output:
        os.path.join(BACT_FASTQ_FIRST_PASS,"{sample}_bacteria_R1.fastq"),
        os.path.join(BACT_FASTQ_FIRST_PASS,"{sample}_bacteria_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem,
        time=180
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
        os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"),
        os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz")
    output:
        os.path.join(VIR_FASTQ_FIRST_PASS,"{sample}_virus_R1.fastq"),
        os.path.join(VIR_FASTQ_FIRST_PASS,"{sample}_virus_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem,
        time=120
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 10239 --include-children  --fastq-output
        """


rule aggr_extraction:
    """aggr"""
    input:
        expand(os.path.join(BACT_FASTQ_FIRST_PASS,"{sample}_bacteria_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(VIR_FASTQ_FIRST_PASS,"{sample}_virus_R1.fastq"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_extraction.txt"),
        time=5
    threads:
        1
    resources:
        mem_mb=SmallJobMem
    shell:
        """
        touch {output[0]}
        """
