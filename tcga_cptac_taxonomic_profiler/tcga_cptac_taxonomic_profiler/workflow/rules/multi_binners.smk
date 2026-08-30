# Additional binners + consensus refinement.
#
# The pipeline previously ran VAMB alone. The v1 MAGs in the manuscript were
# produced by Aviary using CONCOCT (42 bins), MaxBin2 (23) and MetaBAT (16) -
# none of them VAMB - so VAMB-only was a regression on what actually generated
# the published genomes.
#
# All binners here reuse the multi-sample catalogue and per-sample sorted BAMs
# built for VAMB. Multi-sample is the best-performing binning mode in recent
# benchmarks (Nat Commun 2025; Brief Bioinform 2025).
#
# Three algorithm families, so the consensus diversifies the failure mode rather
# than just the tool count:
#   composition + coverage models : MetaBAT2, CONCOCT
#   deep learning embeddings      : VAMB, SemiBin2, COMEBin
#   consensus / refinement        : Binette (scores with CheckM2)

rule metabat2_depth:
    """Contig depth from the sorted BAMs, for MetaBAT2 and CONCOCT."""
    input:
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        depth = os.path.join(METABAT2_RESULTS, 'depth.txt')
    conda:
        os.path.join('..', 'envs', 'metabat2.yaml')
    log:
        os.path.join(LOGS, 'binners', 'metabat2_depth.log')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.med.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.med.cpu
    shell:
        """
        mkdir -p $(dirname {output.depth})
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bams} 2> {log}
        """


rule run_metabat2:
    """MetaBAT2: tetranucleotide frequency + coverage correlation."""
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        depth = os.path.join(METABAT2_RESULTS, 'depth.txt')
    output:
        flag = os.path.join(FLAGS, 'metabat2.flag')
    params:
        outprefix = os.path.join(METABAT2_RESULTS, 'bins', 'bin'),
        min_contig = config.binning.min_contig_length
    conda:
        os.path.join('..', 'envs', 'metabat2.yaml')
    log:
        os.path.join(LOGS, 'binners', 'metabat2.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'metabat2.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p $(dirname {params.outprefix})
        metabat2 --inFile {input.catalogue} --abdFile {input.depth} \
            --outFile {params.outprefix} --minContig {params.min_contig} \
            --numThreads {threads} > {log} 2>&1
        touch {output.flag}
        """


rule run_concoct:
    """CONCOCT: Gaussian mixture model on composition + coverage.

    Produced 42 of the 81 bins in the v1 manuscript - the single largest
    contributor - so it is included for continuity as well as diversity.
    """
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        flag = os.path.join(FLAGS, 'concoct.flag')
    params:
        outdir = CONCOCT_RESULTS,
        min_contig = config.binning.min_contig_length
    conda:
        os.path.join('..', 'envs', 'concoct.yaml')
    log:
        os.path.join(LOGS, 'binners', 'concoct.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'concoct.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p {params.outdir}/bins
        cd {params.outdir}
        zcat {input.catalogue} > catalogue.fna
        cut_up_fasta.py catalogue.fna -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa 2>> {log}
        concoct_coverage_table.py contigs_10K.bed {input.bams} > coverage_table.tsv 2>> {log}
        concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv \
            -b concoct_out --threads {threads} >> {log} 2>&1
        merge_cutup_clustering.py concoct_out_clustering_gt1000.csv > clustering_merged.csv 2>> {log}
        extract_fasta_bins.py catalogue.fna clustering_merged.csv --output_path bins >> {log} 2>&1
        rm -f catalogue.fna contigs_10K.fa
        touch {output.flag}
        """


rule run_semibin2:
    """SemiBin2 multi-sample mode. Among the best performers in recent benchmarks."""
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        flag = os.path.join(FLAGS, 'semibin2.flag')
    params:
        outdir = SEMIBIN2_RESULTS,
        separator = config.binning.separator
    conda:
        os.path.join('..', 'envs', 'semibin2.yaml')
    log:
        os.path.join(LOGS, 'binners', 'semibin2.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'semibin2.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        SemiBin2 multi_easy_bin \
            --input-fasta {input.catalogue} \
            --input-bam {input.bams} \
            --output {params.outdir} \
            --separator {params.separator} \
            --threads {threads} > {log} 2>&1
        touch {output.flag}
        """


rule run_comebin:
    """COMEBin: contrastive multiple-view learning.

    Best-performing deep learning binner in recent benchmarks but also the
    slowest, which is acceptable here since throughput is not the constraint.

    NOTE: run_comebin.sh defaults to -d cuda. Compute nodes are CPU, so the
    device is set explicitly from config; leaving it default would fail.
    """
    input:
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz'),
        bams = expand(os.path.join(VAMB_BAMS, '{sample}_sorted.bam'), sample=SAMPLES)
    output:
        flag = os.path.join(FLAGS, 'comebin.flag')
    params:
        outdir = COMEBIN_RESULTS,
        bamdir = VAMB_BAMS,
        device = config.binning.comebin_device
    conda:
        os.path.join('..', 'envs', 'comebin.yaml')
    log:
        os.path.join(LOGS, 'binners', 'comebin.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'comebin.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p {params.outdir}
        zcat {input.catalogue} > {params.outdir}/catalogue.fna
        run_comebin.sh -a {params.outdir}/catalogue.fna \
            -o {params.outdir} \
            -p {params.bamdir} \
            -t {threads} \
            -d {params.device} > {log} 2>&1
        rm -f {params.outdir}/catalogue.fna
        touch {output.flag}
        """


rule run_binette:
    """Consensus refinement across all five bin sets, scored with CheckM2.

    Binette selects the best-scoring bin per genome across the input bin
    directories, so a diverse candidate set matters more than any single
    binner's performance.
    """
    input:
        vamb = os.path.join(FLAGS, 'vamb.flag'),
        metabat2 = os.path.join(FLAGS, 'metabat2.flag'),
        concoct = os.path.join(FLAGS, 'concoct.flag'),
        semibin2 = os.path.join(FLAGS, 'semibin2.flag'),
        comebin = os.path.join(FLAGS, 'comebin.flag'),
        catalogue = os.path.join(VAMB_CATALOGUE, 'catalogue.fna.gz')
    output:
        flag = os.path.join(FLAGS, 'binette.flag')
    params:
        outdir = BINETTE_RESULTS,
        checkm2_db = config.databases.checkm2,
        bin_dirs = " ".join([
            os.path.join(VAMB_RESULTS, 'bins'),
            os.path.join(METABAT2_RESULTS, 'bins'),
            os.path.join(CONCOCT_RESULTS, 'bins'),
            os.path.join(SEMIBIN2_RESULTS, 'bins'),
            os.path.join(COMEBIN_RESULTS, 'comebin_res', 'comebin_res_bins'),
        ]),
        contamination_weight = config.binning.binette_contamination_weight
    conda:
        os.path.join('..', 'envs', 'binette.yaml')
    log:
        os.path.join(LOGS, 'binners', 'binette.log')
    benchmark:
        os.path.join(BENCHMARKS, 'binners', 'binette.txt')
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time,
        tmpdir = config.tmpdir
    threads:
        config.resources.big.cpu
    shell:
        """
        mkdir -p {params.outdir}
        zcat {input.catalogue} > {params.outdir}/catalogue.fna
        binette --bin_dirs {params.bin_dirs} \
            --contigs {params.outdir}/catalogue.fna \
            --outdir {params.outdir} \
            --checkm2_db {params.checkm2_db} \
            --contamination_weight {params.contamination_weight} \
            --threads {threads} > {log} 2>&1
        rm -f {params.outdir}/catalogue.fna
        touch {output.flag}
        """
