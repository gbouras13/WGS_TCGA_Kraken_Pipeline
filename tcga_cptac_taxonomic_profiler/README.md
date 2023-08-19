# TCGA_CPTAC_Taxonomic_Profiler

Snakemake Pipeline to Mine TCGA and CPTAC WGS data for bacterial reads


Installation
=========

```
git clone "https://github.com/gbouras13/tcga_cptac_taxonomic_profiler.git"
cd tcga_cptac_taxonomic_profiler/
pip install -e .
tcga_cptac_taxonomic_profiler --help
tcga_cptac_taxonomic_profiler install --help
# to run the kraken2/bracken based analysis
tcga_cptac_taxonomic_profiler kraken --help
# to run the mmseqs2 and Uniref50 based analysis
tcga_cptac_taxonomic_profiler mmseqs2 --help
```

Steps:


1. With all BAM files (.bam) files as downloaded from the TCGA/CPTAC/ where-ever you want in the `input_bams` directory, this will extract all unaligned paired end reads (as marked in the BAM) and generate total read counts for each BAM. The unaligned reads will be in the `UNALIGNED_FASTQ` directory.

```
tcga_cptac_taxonomic_profiler extract --input input_bams --output TCGA_output 
```

2. Download the CHM13 `human-t2t-hla.fa` host genome (from [hostile](https://github.com/bede/hostile)) with the following command.

```
tcga_cptac_taxonomic_profiler install_host --database host_genome_db
```

3. Run [trimnami](https://github.com/beardymcjohnface/Trimnami) specifying the directory of FASTQ reads as `--reads`.

Note that you will need to modify the config file so that `--length_required 40` under `qc` and `fastp`. This is because a lot of the TCGA reads are 50 bp.

```
trimnami config
# edit the file so that --length_required 40
trimnami run --reads TCGA_output/UNALIGNED_FASTQ --host host_genome_db/human-t2t-hla.fa fastp --output TCGA_output/trimnami_output --configfile trimnami.config.yaml
```



1. Run the 

You will find the host depleted and trimmed reads in `TCGA_output/trimnami_output/fastp`.

