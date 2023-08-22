"""
Per-sample assemblies for short paired reads
    Take all trimmed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in combine_sample_assemblies.smk

    # taken & modified from hecatomb https://github.com/shandley/hecatomb/blob/main/hecatomb/snakemake/workflow/rules/assembly/paired_end.smk

"""


rule co_assembly:
    """Alternative to cross assembly; assemble everything together in one hit.
    
    Megahit: https://github.com/voutcn/megahit
    """
    input:
        expand(
            os.path.join(INPUT, "{sample}_{r12}.host_rm.fastq.gz"),
            sample=SAMPLES,
            r12=["R1","R2", "S"],
        )
    output:
        assembly = os.path.join(COASSEMBLY_RESULTS, "co_assembly.fasta"),
        graph = os.path.join(COASSEMBLY_RESULTS, "co_assembly_graph.gfa"),
        tmp = temp(os.path.join(COASSEMBLY, "co_assembly_graph.fastg")),
        tar = os.path.join(COASSEMBLY_RESULTS,"coAssembly.tar.zst")
    params:
        r1p = ','.join(expand(os.path.join(INPUT, "{sample}_R1.host_rm.fastq.gz"), sample=SAMPLES)),
        r2p = ','.join(expand(os.path.join(INPUT, "{sample}_R2.host_rm.fastq.gz"), sample=SAMPLES)),
        rs = ','.join(expand(os.path.join(INPUT, "{sample}_S.host_rm.fastq.gz"), sample=SAMPLES)),
        mh_dir = COASSEMBLY,
        mh_int = os.path.join(COASSEMBLY, "intermediate_contigs"),
        params = config.assembly.megahit,
        assembly = os.path.join(COASSEMBLY, "final.contigs.fa"),
        graph = os.path.join(COASSEMBLY, "assembly_graph.gfa"),
    benchmark:
        os.path.join(BENCHMARKS, 'coassembly', "megahit_coassembly.txt")
    log:
        os.path.join(LOGS, 'coassembly', "megahit_coassembly.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join("..", "envs", "megahit.yaml")
    shell:
        """
        if [ -d {params.mh_dir} ]; then
            rm -rf {params.mh_dir}
        fi
        megahit -1 {params.r1p} -2 {params.r2p} -r {params.rs} \
            -o {params.mh_dir} -t {threads} {params.params} &> {log}
        kctg=$(ls -t {params.mh_int}/*.contigs.fa | grep -v final | head -1)
        kmax=$([[ $kctg =~ ([0-9]+) ]] && echo "${{BASH_REMATCH[1]}}")
        megahit_toolkit contig2fastg $kmax $kctg > {output.tmp}
        Bandage reduce {output.tmp} {output.graph}
        cp {params.assembly} {output.assembly}
        tar cf - {params.mh_dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
        rm {log}
        """


rule aggr_coassembly:
    input:
        assembly = os.path.join(COASSEMBLY_RESULTS, "co_assembly.fasta"),
        graph = os.path.join(COASSEMBLY_RESULTS, "co_assembly_graph.gfa"),
        tar = os.path.join(COASSEMBLY_RESULTS,"coAssembly.tar.zst")
    output:
        os.path.join(FLAGS, "aggr_coassembly.flag")
    threads:
        1
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    shell:
        """
        touch {output[0]}
        """
