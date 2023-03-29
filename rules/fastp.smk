        
rule concat__pre_fastp:
    """concat before fastp"""
    input:
        os.path.join(BACT_FASTQ_FIRST_PASS,"{sample}_bacteria_R1.fastq"),
        os.path.join(BACT_FASTQ_FIRST_PASS,"{sample}_bacteria_R2.fastq"),
        os.path.join(VIR_FASTQ_FIRST_PASS,"{sample}_virus_R1.fastq"),
        os.path.join(VIR_FASTQ_FIRST_PASS,"{sample}_virus_R2.fastq"),
        os.path.join(FUNGI_FASTQ_FIRST_PASS,"{sample}_fungi_R1.fastq"),
        os.path.join(FUNGI_FASTQ_FIRST_PASS,"{sample}_fungi_R2.fastq")
    output:
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_R1.fastq"),
        os.path.join(CONCAT_FASTQ,"{sample}_bacteria_virus_fungi_R2.fastq")
    resources:
        mem_mb=SmallJobMem,
        time=5
    threads:
        1
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
        mem_mb=16000,
        time=2800
    threads:
        16
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
        --thread {threads}
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
        time=5
    threads:
        1
    shell:
        """
        touch {output[0]}
        """


rule create_host_index:
    """Step 02. Create the minimap2 index for mapping to the host; this will save time."""
    input:
        HOSTFA,
    output:
        HOSTINDEX
    benchmark:
        os.path.join(BENCH, "create_host_index.txt")
    log:
        os.path.join(STDERR, 'create_host_index.log')
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -t {threads} -d {output} <(cat {input}) 2> {log}
        rm {log}
        """

rule host_removal_mapping:
    """Preprocessing step 02a: Host removal: mapping to host.
    
    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = os.path.join(TMPDIR, "p01", "{sample}_R1.s1.out.fastq"),
        r2 = os.path.join(TMPDIR, "p01", "{sample}_R2.s1.out.fastq"),
        host = HOSTINDEX
    output:
        r1 = temp(os.path.join(TMPDIR, "p02", "{sample}_R1.unmapped.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p02", "{sample}_R2.unmapped.fastq")),
        s = temp(os.path.join(TMPDIR, "p02", "{sample}_R1.unmapped.singletons.fastq")),
        o = temp(os.path.join(TMPDIR, "p02", "{sample}_R1.other.singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "host_removal_mapping.{sample}.txt")
    log:
        mm = os.path.join(STDERR, "host_removal_mapping.{sample}.minimap.log"),
        sv = os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsView.log"),
        fq = os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} --secondary=no {input.host} {input.r1} {input.r2} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fastq -NO -1 {output.r1} -2 {output.r2} -0 {output.o} -s {output.s} 2> {log.fq}
        rm {log.mm} {log.sv} {log.fq}
        """