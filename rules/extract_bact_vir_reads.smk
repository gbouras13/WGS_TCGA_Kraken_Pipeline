rule extract_bact_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R1.unmapped.fastq"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R2.unmapped.fastq")
    output:
        os.path.join(BACT_FASTQ,"{sample}_bacteria_R1.fastq"),
        os.path.join(BACT_FASTQ,"{sample}_bacteria_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=180
    threads:
        1
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 2 --include-children  --fastq-output
        """

rule extract_virus_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R1.unmapped.fastq"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R2.unmapped.fastq")
    output:
        os.path.join(VIR_FASTQ,"{sample}_virus_R1.fastq"),
        os.path.join(VIR_FASTQ,"{sample}_virus_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=120
    threads:
        1
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 10239 --include-children  --fastq-output
        """

        
rule extract_hpv_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R1.unmapped.fastq"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R2.unmapped.fastq")
    output:
        os.path.join(VIR_FASTQ,"{sample}_hpv_R1.fastq"),
        os.path.join(VIR_FASTQ,"{sample}_hpv_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=120
    threads:
        1
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 151340 --include-children  --fastq-output
        """
        
rule extract_fung_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN,"{sample}.kraken.txt"),
        os.path.join(KRAKEN,"{sample}.kraken.rep"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R1.unmapped.fastq"),
        os.path.join(HOST_UNMAPPED_FASTQ, "{sample}_R2.unmapped.fastq")
    output:
        os.path.join(FUNGI_FASTQ,"{sample}_fungi_R1.fastq"),
        os.path.join(FUNGI_FASTQ,"{sample}_fungi_R2.fastq")
    params:
        os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
    conda:
        os.path.join('..', 'envs','kraken2.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=180
    threads:
        1
    shell:
        """
        python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 4751 --include-children  --fastq-output
        """
        

rule aggr_extraction:
    """aggr"""
    input:
        expand(os.path.join(BACT_FASTQ,"{sample}_bacteria_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(FUNGI_FASTQ,"{sample}_fungi_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(VIR_FASTQ,"{sample}_virus_R1.fastq"), sample = SAMPLES),
        expand(os.path.join(VIR_FASTQ,"{sample}_hpv_R1.fastq"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_extraction.txt")
    resources:
        mem_mb=SmallJobMem,
        time=5
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
