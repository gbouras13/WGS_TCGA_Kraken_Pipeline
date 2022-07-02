rule fastp:
    """use fastp to qc files"""
    input:
        os.path.join(TMP,"{sample}_bacteria_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_bacteria_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_bacteria_fastp_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_bacteria_fastp_R2.fastq.gz")
    log:
        os.path.join(LOGS,"{sample}.fastp.log")
    conda:
        os.path.join('..', 'envs','fastp.yaml')
    params:
        os.path.join(CONTAMINANTS, 'vector_contaminants.fa')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
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
        --dedup --dup_calc_accuracy 4   --trim_poly_x 
        --thread {threads}
        """
        
rule aggr_fastp:
    """aggr"""
    input:
        expand(os.path.join(TMP,"{sample}_bacteria_fastp_R1.fastq.gz"), sample = SAMPLES),
    output:
        os.path.join(LOGS, "aggr_fastp.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem
    shell:
        """
        touch {output[0]}
        """
        