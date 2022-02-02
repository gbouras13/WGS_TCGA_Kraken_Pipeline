# WGS_TCGA_Pipeline
Snakemake Pipeline to Mine TCGA WGS Data for Contaminant Reads

#### This is a Work in Progress.

* Some code (namely, samples.smk) and inspiration for the general structure has been borrowed from https://github.com/shandley/hecatomb
* Inputs required are the relevant Bam files - ideally in the Bams/ directory, but can be anywhere.
* Only requirement is that snakemake be in the $PATH.

# Usage

1. Download the Kraken DB
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

* with Slurm:

```console
snakemake -c 16 -s wgs_runner.smk --use-conda --config Reads=Bams/
```
