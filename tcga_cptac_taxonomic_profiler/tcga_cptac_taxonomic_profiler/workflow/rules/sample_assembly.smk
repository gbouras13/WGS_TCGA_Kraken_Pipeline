"""
Per-sample assemblies for short paired reads
    # taken & modified from hecatomb https://github.com/shandley/hecatomb/blob/main/hecatomb/snakemake/workflow/rules/assembly/paired_end.smk
    # and aviary https://github.com/rhysnewell/aviary/blob/master/aviary/modules/assembly/assembly.smk
"""

rule individual_sample_assembly:
    """Individual sample assemblies
        Use metaspades as Vijini says it is better than megahit :)
    """
    input:
        r1p = os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz"), 
        r2p = os.path.join(INPUT, "{sample}_R2.host_rm.fastq.gz"), 
        rs = os.path.join(INPUT, "{sample}_S.host_rm.fastq.gz")
    output:
        fasta = os.path.join(SAMPLE_ASSEMBLIES, '{sample}', 'scaffolds.fasta') 
    params:
        max_memory = config.resources.big.mem,
        assembly_dir = os.path.join(SAMPLE_ASSEMBLIES, '{sample}'),
        tmpdir = TMPDIR
    benchmark:
        os.path.join(BENCHMARKS, 'sample_assembly', "{sample}_spades_assembly.txt")
    log:
        os.path.join(LOGS, 'sample_assembly', "{sample}_spades_assembly.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join("..", "envs", "spades.yaml")
    shell:
        """
        # aviary

        minimumsize=200000;
        actualsize=$(stat -c%s {input.r1p});

        if [ $actualsize -ge $minimumsize ]
        then
            spades.py --restart-from last --checkpoints all --memory {params.max_memory} --meta -1 {input.r1p} -2 {input.r2p} -s {input.rs}  \
            -o {params.assembly_dir} -t {threads}  -k auto --tmp-dir {params.tmpdir} 2> {log} 
        else
            touch {output.fasta}

        """


rule aggr_sample_assembly:
    input:
        expand(os.path.join(SAMPLE_ASSEMBLIES, '{sample}', 'scaffolds.fasta'), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_sample_assembly.flag")
    threads:
        1
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """
