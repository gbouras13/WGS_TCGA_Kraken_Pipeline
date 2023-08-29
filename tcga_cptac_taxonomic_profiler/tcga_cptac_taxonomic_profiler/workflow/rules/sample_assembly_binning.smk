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
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
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
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
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
        minfasta = config.binning.minfasta,
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

"""
checkm2
"""

# Evaluate in which samples bins were reconstructed
checkpoint samples_with_bins:
    input:        
        outtouch = os.path.join(FLAGS, 'vamb.flag')
    output:
        sammples_with_bins = os.path.join(VAMB_CATALOGUE, 'samples_with_bins.txt')
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    log:
        os.path.join(LOGS, 'vamb', "samples_with_bins.log")
    params:
        outdir = VAMB_RESULTS
    threads:
        1
    shell:
        """
        find {params.outdir}/bins/*/ -type d ! -empty |sed 's=.*bins/==g'  |sed 's=/==g'  > {output.sammples_with_bins}
        """

def samples_with_bins_f(wildcards):
    # decision based on content of output file
    with checkpoints.samples_with_bins.get().output[0].open() as f:
        samples_with_bins = [sample.strip() for sample in f.readlines()]
        samples_with_bins_paths=expand(os.path.join(CHECKM2_RESULTS,"flags/checkm2_all_{sample}_bins_finished.flag"),sample=samples_with_bins)
        return samples_with_bins_paths


# Run CheckM2 for each sample with bins        
rule run_checkm2_per_sample_all_bins:
    output:
        outtouch=os.path.join(CHECKM2_RESULTS,"flags/checkm2_all_{sample}_bins_finished.flag")
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.big.cpu
    params:
        vamb_dir = VAMB_RESULTS,
        checkm2_dir = CHECKM2_RESULTS,
        database = config.databases.checkm2
    log:
        os.path.join(LOGS, 'vamb', "{sample}_checkm2.log")
    conda: 
        os.path.join("..", "envs", "checkm2.yaml")
    shell:
        """
        checkm2 predict --threads {threads} --input {params.vamb_dir}/bins/{wildcards.sample}/*.fna --database_path {params.database} --output-directory {params.checkm2_dir}/{wildcards.sample} --force
        touch {output.outtouch}
        """


def samples_with_bins(wildcards):
    # decision based on content of output file
    with checkpoints.samples_with_bins.get().output[0].open() as f:
        samples_with_bins = [sample.strip() for sample in f.readlines()]
        return samples_with_bins


# this rule will be executed when all CheckM2 runs per sample finish
rule cat_checkm2_all:
    input:
        samples_with_bins_f
    output: 
        outtouch = os.path.join(FLAGS, 'checkm2.flag')
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    threads:
        1
    shell:
        "touch {output}"



rule get_hq_bins:
    """
    get_hq_bins
    """
    input:
        samples_with_bins_f,
        outtouch = os.path.join(FLAGS, 'checkm2.flag')
    params:
        samples = samples_with_bins,
        checkm2_directory = CHECKM2_RESULTS,
        vamb_bin_dir = VAMB_RESULTS,
        combined_mag_directory = ALL_MAGS
    output:
        checm2_combo = os.path.join(CHECKM2_RESULTS, "combined_check2_quality_report_hq.tsv")
    benchmark:
        os.path.join(BENCHMARKS, 'vamb', "get_hq_bins.txt")
    log:
        os.path.join(LOGS, 'vamb', "get_hq_bins.log")
    resources:
        mem_mb = config.resources.sml.mem,
        time = config.resources.sml.time
    threads:
        config.resources.sml.cpu
    conda:
        os.path.join("..", "envs", "biopython.yaml")
    script:
        '../scripts/get_med_hq_bins.py'



"""
gtdbtk
# """

rule gtdbtk_ani:
    """
    generate mash db
    """
    input:
        samples_with_bins_f,
        outtouch = os.path.join(FLAGS, 'checkm2.flag'),
    output:
        out_tsv = os.path.join(GTDB_MASH_OUTDIR, 'gtdbtk.ani_summary.tsv') 
    threads:
        config.resources.big.cpu
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.med.time
    conda:
        os.path.join("..", "envs", "gtdbtk.yaml")
    benchmark:
        os.path.join(BENCHMARKS, 'gtdb', "gtdbtk_ani.txt")
    log:
        os.path.join(LOGS, 'gtdb', "gtdbtk_ani.log")
    params:
        combined_mag_directory = ALL_MAGS,
        mashdir=GTDB_MASH_OUTDIR,
        GTDBTK_DATA_PATH = config.databases.gtdb
    shell:
        """
        export GTDBTK_DATA_PATH={params.GTDBTK_DATA_PATH}

        gtdbtk ani_rep --genome_dir {params.combined_mag_directory} --out_dir {params.mashdir} --cpus {threads}

        """


rule gtdbtk_classify_wf:
    input:
        samples_with_bins_f,
        outtouch = os.path.join(FLAGS, 'checkm2.flag'),
        out_tsv = os.path.join(GTDB_MASH_OUTDIR, 'gtdbtk.ani_summary.tsv') 
    output:
        outtouch = os.path.join(FLAGS, 'gtdb.flag')
    threads:
        config.resources.big.cpu
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.med.time
    conda:
        os.path.join("..", "envs", "gtdbtk.yaml")
    benchmark:
        os.path.join(BENCHMARKS, 'gtdb', "gtdbtk_classify_wf.txt")
    log:
        os.path.join(LOGS, 'gtdb', "gtdbtk_classify_wf.log")
    params:
        combined_mag_directory = ALL_MAGS,
        outdir=GTDB_OUTDIR,
        GTDBTK_DATA_PATH = config.databases.gtdb,
        mashdir=GTDB_MASH_OUTDIR,
    shell:
        """
        export GTDBTK_DATA_PATH={params.GTDBTK_DATA_PATH}

        gtdbtk classify_wf --genome_dir {params.combined_mag_directory}  --out_dir {params.outdir} --cpus {threads}   --mash_db {params.mashdir}  --force

        touch {output.outtouch}
        """


# # Evaluate in which mags were created
# checkpoint get_mags:
#     input:        
#         outtouch = os.path.join(FLAGS, 'mag.flag')
#     output:
#         mags = os.path.join(BAKTA, 'mag_names.txt')
#     resources:
#         mem_mb = config.resources.sml.mem,
#         time = config.resources.sml.time
#     log:
#         os.path.join(LOGS,  "mag_names.log")
#     params:
#         outdir = ALL_MAGS
#     threads:
#         1
#     shell:
#         """
#         find {params.outdir}/* -type d ! -empty |sed 's=.*bins/==g'  |sed 's=/==g'  > {output.mags}
#         """

# def all_mags(wildcards):
#     # decision based on content of output file
#     with checkpoints.all_mags.get().output[0].open() as f:
#         all_mags = [sample.strip() for sample in f.readlines()]
#         return all_mags

# rule run_bakta:
#     """
#     generate mash db
#     """
#     input:
#         all_mags,
#         mag =  os.path.join(ALL_MAGS, '{mag}.fna')
#     output:
#         out_tsv = os.path.join(BAKTA, '{mag}', '{mag}.tsv') 
#     threads:
#         config.resources.med.cpu
#     resources:
#         mem_mb = config.resources.med.mem,
#         time = config.resources.med.time
#     conda:
#         os.path.join("..", "envs", "bakta.yaml")
#     benchmark:
#         os.path.join(BENCHMARKS, 'bakta', "bakta_{mag}.txt")
#     log:
#         os.path.join(LOGS, 'bakta', "bakta_{mag}.log")
#     params:
#         outdir = os.path.join(BAKTA, '{mag}')
#         db=config.databases.bakta,
#     shell:
#         """

#         bakta --db {params.db} --output {params.outdir} -f  {input.mag} 

#         """

# rule aggr_bakta:
#     """
#     generate mash db
#     """
#     input:
#         all_mags,
#         tsvs = expand(os.path.join(BAKTA, '{mag}', '{mag}.tsv') , mag=all_mags),
#         outtouch = os.path.join(FLAGS, 'checkm2.flag')
#     output:
#         outtouch = os.path.join(FLAGS, 'bakta.flag'),
#     threads:
#         config.resources.sml.cpu
#     resources:
#         mem_mb = config.resources.sml.mem,
#         time = config.resources.sml.time
#     shell:
#         """
#         touch {output.outtouch}

#         """
