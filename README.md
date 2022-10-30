# WGS_TCGA_Pipeline
Snakemake Pipeline to Mine WGS Data for Contaminant Reads

**This is a Work in Progress.**

* This pipeline is not actually TCGA specific, but it does require as input WGS short reads reads mapped against a host genome in .bam format - it should actually work for all non bacterial/virus hosts.
* Some code (namely, the sample parsing in samples.smk and decontamination in fastp.smk) has been borrowed and modified from https://github.com/shandley/hecatomb.
* The only input required are the relevant Bam files, which all must be placed in a certain directory specified as "Bams"
* Only software requirement is conda, and that snakemake be in the $PATH. The rest of the required programs should install via conda
* extract_reads.py has been taken from KrakenTools https://github.com/jenniferlu717/KrakenTools#extract_kraken_readspy

# Usage

1. Download the Kraken2 DB
* This needs to only be run once

```console
snakemake -c 1 -s DownloadDB.smk
```

2. Next the unaligned reads need to be extracted from the bams.

* This was split from the main pipeline due to the massive data size of the input data files (and so can be run independently of the main pipeline, allowing you to delete the raw BAMs once the unaligned reads have been extracted).
* All you need to specify is the Bams directory, and an output directory.

```console
snakemake -c 1 -s extract_unaligned_fastq --use-conda --config Bams=Bams/ Output=TCGA_Output/ 
```

3. Run the pipeline

```console
snakemake -c 16 -s wgs_runner.smk --use-conda --config Output=my_output_dir/
```

Other Notes
======

* For offline use (e.g. Adelaide Uni HPC) - the conda envs need to be installed first on the login node in the pipeline directory

```console
snakemake -c 1 -s wgs_runner.smk --use-conda --config Output=my_output_dir/ --conda-create-envs-only --conda-frontend conda
```

* It it highly recommended that you run this piepline With a Slurm profile (see https://snakemake.readthedocs.io/en/stable/executing/cli.html https://github.com/Snakemake-Profiles/slurm https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated)
