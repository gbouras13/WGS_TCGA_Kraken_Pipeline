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

rule fastp_trim:
    """remove adapters etc """
    input:
        os.path.join(TMP,"{sample}_bacteria_R1.fastq"),
        os.path.join(TMP,"{sample}_bacteria_R2.fastq")
    output:
        os.path.join(TMP,"{sample}_bacteria_trim_R1.fastq"),
        os.path.join(TMP,"{sample}_bacteria_trim_R2.fastq")
    log:
        os.path.join(LOGS,"{sample}_fastp.log")
    conda:
        os.path.join('..', 'envs','assembly.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -w {threads}
        """


rule megahit:
    """run megahit."""
    input:
        os.path.join(TMP,"{sample}_bacteria_trim_R1.fastq"),
        os.path.join(TMP,"{sample}_bacteria_trim_R2.fastq")
    output:
        os.path.join(MEGAHIT, "{sample}", "final.contigs.fa")
    params:
        directory(os.path.join(MEGAHIT, '{sample}'))
    log:
        os.path.join(LOGS,"{sample}.megahit.log")
    conda:
        os.path.join('..', 'envs','assembly.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
	    rm -rf {params[0]}
        megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
        """

# rule extract_fuso_fastqs:
#     """Extract Fuso Fastas - taxid 848."""
#     input:
#         os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
#         os.path.join(TMP,"{sample}_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_R2.fastq.gz")
#     output:
#         os.path.join(TMP,"{sample}_fuso_R1.fastq"),
#         os.path.join(TMP,"{sample}_fuso_R2.fastq")
#     params:
#         os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
#     log:
#         os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
#         -o {output[0]} -o2 {output[1]} -r {input[1]} -t 848 --include-children  --fastq-output
#         """
        
# rule extract_strep_fastqs:
#     """Extract Strep Fastas - taxid 1301."""
#     input:
#         os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
#         os.path.join(TMP,"{sample}_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_R2.fastq.gz")
#     output:
#         os.path.join(TMP,"{sample}_strep_R1.fastq"),
#         os.path.join(TMP,"{sample}_strep_R2.fastq")
#     params:
#         os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
#     log:
#         os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
#         -o {output[0]} -o2 {output[1]} -r {input[1]} -t 1301 --include-children  --fastq-output
#         """

# rule extract_prevotella_fastqs:
#     """Extract Fuso Fastas - taxid 838."""
#     input:
#         os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
#         os.path.join(TMP,"{sample}_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_R2.fastq.gz")
#     output:
#         os.path.join(TMP,"{sample}_prev_R1.fastq"),
#         os.path.join(TMP,"{sample}_prev_R2.fastq")
#     params:
#         os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
#     log:
#         os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
#         -o {output[0]} -o2 {output[1]} -r {input[1]} -t 838 --include-children  --fastq-output
#         """

# rule extract_citrobacter_fastqs:
#     """Extract citro Fastas - taxid 544."""
#     input:
#         os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
#         os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
#         os.path.join(TMP,"{sample}_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_R2.fastq.gz")
#     output:
#         os.path.join(TMP,"{sample}_citro_R1.fastq"),
#         os.path.join(TMP,"{sample}_citro_R2.fastq")
#     params:
#         os.path.join(KRAKENTOOLSDIR, 'extract_kraken_reads.py')
#     log:
#         os.path.join(LOGS,"{sample}.extract_kraken_reads.log")
#     conda:
#         os.path.join('..', 'envs','kraken2.yaml')
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         python3 {params[0]} -k {input[0]} -s1 {input[2]} -s2 {input[3]} \
#         -o {output[0]} -o2 {output[1]} -r {input[1]} -t 544 --include-children  --fastq-output
#         """

# rule megahit_citro:
#     """run megahit_fuso."""
#     input:
#         os.path.join(TMP,"{sample}_citro_R1.fastq"),
#         os.path.join(TMP,"{sample}_citro_R2.fastq")
#     output:
#         os.path.join(MEGAHIT, "{sample}_Citro", "final.contigs.fa")
#     params:
#         directory(os.path.join(MEGAHIT, '{sample}_Citro'))
#     log:
#         os.path.join(LOGS,"{sample}.megahit_citro.log")
#     conda:
#         os.path.join('..', 'envs','assembly.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
# 	    rm -rf {params[0]}
#         megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
#         """
	
# rule megahit_fuso:
#     """run megahit_fuso."""
#     input:
#         os.path.join(TMP,"{sample}_fuso_R1.fastq"),
#         os.path.join(TMP,"{sample}_fuso_R2.fastq")
#     output:
#         os.path.join(MEGAHIT, "{sample}_Fuso", "final.contigs.fa")
#     params:
#         directory(os.path.join(MEGAHIT, '{sample}_Fuso'))
#     log:
#         os.path.join(LOGS,"{sample}.megahit_fuso.log")
#     conda:
#         os.path.join('..', 'envs','assembly.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
# 	    rm -rf {params[0]}
#         megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
#         """

# rule megahit_prev:
#     """run megahit_prev."""
#     input:
#         os.path.join(TMP,"{sample}_prev_R1.fastq"),
#         os.path.join(TMP,"{sample}_prev_R2.fastq")
#     output:
#         os.path.join(MEGAHIT, "{sample}_Prev", "final.contigs.fa")
#     params:
#         directory(os.path.join(MEGAHIT, '{sample}_Prev'))
#     log:
#         os.path.join(LOGS,"{sample}.megahit_prev.log")
#     conda:
#         os.path.join('..', 'envs','assembly.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
# 	    rm -rf {params[0]}
#         megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
#         """

# rule megahit_strep:
#     """run megahit_strep."""
#     input:
#         os.path.join(TMP,"{sample}_strep_R1.fastq"),
#         os.path.join(TMP,"{sample}_strep_R2.fastq")
#     output:
#         os.path.join(MEGAHIT, "{sample}_Strep", "final.contigs.fa")
#     params:
#         directory(os.path.join(MEGAHIT, '{sample}_Strep'))
#     log:
#         os.path.join(LOGS,"{sample}.megahit_strep.log")
#     conda:
#         os.path.join('..', 'envs','assembly.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         rm -rf {params[0]}
#         megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
#         """


rule extract_clust_1_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_clust_1_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_1_R2.fastq")
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
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 39948 836 194 1378 86331 543311 848 157 838 1314 195950 32067 1016 970 1164882 816 44259 724 1654 482 416916 538 --include-children  --fastq-output
        """

rule fastp_trim_clust_1:
    """remove adapters etc """
    input:
        os.path.join(TMP,"{sample}_clust_1_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_1_R2.fastq")
    output:
        os.path.join(TMP,"{sample}_clust_1_trim_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_1_trim_R2.fastq")
    log:
        os.path.join(LOGS,"{sample}_fastp.log")
    conda:
        os.path.join('..', 'envs','assembly.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -w {threads}
        """

rule megahit_clust_1:
    """run megahit."""
    input:
        os.path.join(TMP,"{sample}_clust_1_trim_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_1_trim_R2.fastq")
    output:
        os.path.join(MEGAHIT, "{sample}_clust_1", "final.contigs.fa")
    params:
        directory(os.path.join(MEGAHIT, '{sample}_clust_1'))
    log:
        os.path.join(LOGS,"{sample}.megahit.log")
    conda:
        os.path.join('..', 'envs','assembly.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
	    rm -rf {params[0]}
        megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
        """

rule extract_clust_2_fastqs:
    """Extract Fastas."""
    input:
        os.path.join(KRAKEN_S,"{sample}.kraken.txt"),
        os.path.join(KRAKEN_S,"{sample}.kraken.rep"),
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_clust_2_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_2_R2.fastq")
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
        -o {output[0]} -o2 {output[1]} -r {input[1]} -t 544 1822464 1910952 1910951 125216 53335 1827 10509 620 234 10662 32008 2701 1883 1280 10294 165696 10699 2759736 1912216 48736 1763 13687 2742598 1350 286 333750 590 162289 32257 --include-children  --fastq-output
        """

rule fastp_trim_clust_2:
    """remove adapters etc """
    input:
        os.path.join(TMP,"{sample}_clust_2_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_2_R2.fastq")
    output:
        os.path.join(TMP,"{sample}_clust_2_trim_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_2_trim_R2.fastq")
    log:
        os.path.join(LOGS,"{sample}_fastp.log")
    conda:
        os.path.join('..', 'envs','assembly.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -w {threads}
        """

rule megahit_clust_2:
    """run megahit."""
    input:
        os.path.join(TMP,"{sample}_clust_2_trim_R1.fastq"),
        os.path.join(TMP,"{sample}_clust_2_trim_R2.fastq")
    output:
        os.path.join(MEGAHIT, "{sample}_clust_2", "final.contigs.fa")
    params:
        directory(os.path.join(MEGAHIT, '{sample}_clust_2'))
    log:
        os.path.join(LOGS,"{sample}.megahit.log")
    conda:
        os.path.join('..', 'envs','assembly.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
	    rm -rf {params[0]}
        megahit -1 {input[0]} -2 {input[1]} -o {params[0]}
        """


rule aggr_assembly:
    """aggr"""
    input:
        expand(os.path.join(MEGAHIT, "{sample}", "final.contigs.fa"), sample = SAMPLES),
        expand(os.path.join(MEGAHIT, "{sample}_clust_1", "final.contigs.fa"), sample = SAMPLES),
        expand(os.path.join(MEGAHIT, "{sample}_clust_2", "final.contigs.fa"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_assembly.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
