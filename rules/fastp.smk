        
rule concat__pre_fastp:
    """concat before fastp"""
    input:
        os.path.join(BACT_FASTQ_FIRST_PASS,"{sample}_bacteria_R1.fastq.gz"),
        os.path.join(BACT_FASTQ_FIRST_PASS,"{sample}_bacteria_R2.fastq.gz"),
        os.path.join(VIR_FASTQ_FIRST_PASS,"{sample}_virus_R1.fastq.gz"),
        os.path.join(VIR_FASTQ_FIRST_PASS,"{sample}_virus_R2.fastq.gz"),
        os.path.join(FUNGI_FASTQ_FIRST_PASS,"{sample}_fungi_R1.fastq.gz"),
        os.path.join(FUNGI_FASTQ_FIRST_PASS,"{sample}_fungi_R2.fastq.gz")
    output:
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_R1.fastq.gz"),
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_R2.fastq.gz")
    resources:
        mem_mb=SmallJobMem,
        time=5,
        th=1
    shell:
        """
        cat {input[0]} {input[2]} {input[4]} > {output[0]}
        cat {input[1]} {input[3]} {input[4]} > {output[1]}
        """

rule fastp_all:
    """use fastp to qc files"""
    input:
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_R1.fastq"),
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_R2.fastq")
    output:
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_fastp_R1.fastq.gz"),
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_fastp_R2.fastq.gz")
    conda:
        os.path.join('..', 'envs','fastp.yaml')
    params:
        os.path.join(CONTAMINANTS, 'vector_contaminants.fa')
    resources:
        mem_mb=BigJobMem,
        th=32
    shell:
        """
        fastp -i {input[0]} \
        -I {input[1]} \
        -o {output[0]} \
        -O {output[1]}  \
        -z 1 \
        --length_required 40 \
        --adapter_fasta {params[0]} \
        --cut_tail --cut_tail_window_size 25 --cut_tail_mean_quality 15  \
        --dedup --dup_calc_accuracy 4   --trim_poly_x \
        --thread {resources.th}
        """
        
rule aggr_fastp:
    """aggr"""
    input:
        expand(os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_fastp_R1.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_fastp_R2.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_fastp.txt")
    resources:
        mem_mb=SmallJobMem,
        time=5,
        th=1
    shell:
        """
        touch {output[0]}
        """
        
