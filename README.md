# WGS_TCGA_Pipeline
Snakemake Pipeline to Mine TCGA WGS Data for Contaminant Reads

#### This is a Work in Progress.

* Some code (namely, samples.smk) and inspiration for the general structure has been borrowed from https://github.com/shandley/hecatomb
* Inputs required are the relevant Bam files - ideally in the Bams/ directory, but can be anywhere (the directory they are in must be specified with Reads={directory}).
* Only requirement is that snakemake be in the $PATH.

# Usage

1. Download the Kraken2 DB
* This needs to only be run once

```console
snakemake -c 1 -s DownloadDB.smk
```

2. For offline use (e.g. Adelaide Uni HPC) - the conda envs need to be installed first from the login node

```console
snakemake -c 1 -s wgs_runner.smk --use-conda --config Reads=Bams/ --conda-create-envs-only --conda-frontend conda
```

3. Run the pipeline

```console
snakemake -c 16 -s wgs_runner.smk --use-conda --config Reads=Bams/
```

* with a Slurm profile (see https://snakemake.readthedocs.io/en/stable/executing/cli.html https://github.com/Snakemake-Profiles/slurm https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated)
* You will need to cd to the pipeline directory in your jobscript before running if you want to run this offline (to use the premade conda envs)

```console
snakemake -c 16 -s wgs_runner.smk --config Reads=Bams/  --profile wgs_tcga
```
