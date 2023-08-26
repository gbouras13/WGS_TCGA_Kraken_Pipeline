"""
sample assembly binning will use VAMB
https://github.com/RasmussenLab/vamb
"""

rule concatenate_sample_assemblies:
    """
    concatenate assemblies
    """

    input:
        fastas = expand(os.path.join(SAMPLE_ASSEMBLIES, '{sample}', 'contigs.fasta'), sample=SAMPLES)
    params:
        min_contig_length = config.binning.min_contig_length,
        samples = SAMPLES,
        separator = ":"
    output:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        csv_path = os.path.join(VAMB_CATALOGUE, 'sample_bin_name.csv')
    benchmark:
        os.path.join(BENCHMARKS, 'vamb', "concatenate_sample_assemblies.txt")
    log:
        os.path.join(LOGS, 'vamb', "concatenate_sample_assemblies.log")
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join("..", "envs", "biopython.yaml")
    script:
        '../scripts/concatenate_save_samples.py'


rule index_catalogue:
    """
    index catalogue
    """
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz')
    output:
        index = os.path.join(VAMB_CATALOGUE, 'catalogue.mmi')
    benchmark:
        os.path.join(BENCHMARKS, 'vamb', "index_catalogue.txt")
    log:
        os.path.join(LOGS, 'vamb', "index_catalogue.log")
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join("..", "envs", "minimap2.yaml")
    shell:
        """
        minimap2 -d {output.index} {input.catalogue}; # make index
        """

rule read_mapping:
    """
    maps reads to the contig catalogue
    """
    input:
        r1p = os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz"), 
        r2p = os.path.join(INPUT, "{sample}_R2.host_rm.fastq.gz"), 
        index = os.path.join(VAMB_CATALOGUE, 'catalogue.mmi')
    output:
        bam = os.path.join(VAMB_BAMS, '{sample}.bam')
    benchmark:
        os.path.join(BENCHMARKS, 'vamb', "{sample}_read_mapping.txt")
    log:
        os.path.join(LOGS, 'vamb', "{sample}_read_mapping.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join("..", "envs", "minimap2.yaml")
    shell:
        """

        minimap2 -t {threads} -ax sr {input.index}  {input.r1p} {input.r2p}  | 
        samtools view -F 3584 -b --threads {threads} > {output.bam}

        """

rule vamb_bam_Sprt:
    input:
        bam = os.path.join(VAMB_BAMS, '{sample}.bam')
    output:
        bam = os.path.join(VAMB_BAMS, '{sample}_sorted.bam') 
    benchmark:
        os.path.join(BENCHMARKS, 'vamb', "{sample}_read_mapping_sort.txt")
    log:
        os.path.join(LOGS, 'vamb', "{sample}_read_mapping_sort.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join("..", "envs", "minimap2.yaml")
    shell:
        """
        samtools sort {input.bam} --threads {threads} -o {output.bam} 2> {log}
        """


rule run_vamb:
    """
    maps reads to the contig catalogue
    needs vamb 4.1.3 in Path 
    """
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        outtouch = os.path.join(FLAGS, 'vamb.flag')
    benchmark:
        os.path.join(BENCHMARKS, 'vamb', "run_vamb.txt")
    log:
        os.path.join(LOGS, 'vamb', "run_vamb.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    params:
        bams = ' '.join(expand(os.path.join(VAMB_BAMS, '{sample}.bam'), sample=SAMPLES)),
        outdir = VAMB_RESULTS,
        separator = config.binning.separator,
        minfasta = config.binning.minfasta
        min_contig_length = config.binning.min_contig_length
    threads:
        config.resources.big.cpu
    shell:
        # you may have trouble running this - best not to use a profile for this rule
        # By default, Vamb does not output any FASTA files of the bins. 
        # In the examples below, the option --minfasta 200000 is set, meaning that all bins with a size of 200 kbp or more will be output as FASTA files.
        """

        vamb --outdir {params.outdir} --fasta {input.catalogue} --bamfiles {input.bams} -o {params.separator} -m {params.min_contig_length} --minfasta {params.minfasta}
        touch {output.outtouch}
        
        """




