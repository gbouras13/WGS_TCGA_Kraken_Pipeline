# TCGA_CPTAC_Taxonomic_Profiler

Snakemake Pipeline to Mine TCGA and CPTAC HNSC (or any other really) Human WGS data for non-host reads

# Installation

* I recommend installing into a virtual environment (e.g. one that is made by conda/mamba)

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

# Database Downloads

*  You will need to download the relevant databases and change the paths in the config file to where they live on your system:

* You can find these in `tcga_cptac_taxonomic_profiler/config/config.yaml`
* You won't need to do the host database as it is done in step 2.

```
databases:
    host_db: '~/host_genome_db'
    kraken: '~/kraken2/k2_pluspf_20230605'
    mmseqs2:
        uniref50: '~/metagenomics_dbs/mmseqs2/UniRef50/UniRef50'
        uniref90: '~/metagenomics_dbs/mmseqs2/UniRef90/UniRef90'
    checkm2: '/~/metagenomics_dbs/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd'
    gtdb: '~/GTDB/release214/'
    bakta: '~/metagenomics_dbs/bakta/db'
```

* **Kraken 2**

```
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230605.tar.gz
tar -xzf k2_pluspf_20230605.tar.gz
```

* I used the 5 June 2023 PlusPF DB when I did the analysis, there are newer ones now
* You can find plenty of options here https://benlangmead.github.io/aws-indexes/k2
* **As long as you use one with the human genome in it, I am sure the results will be similar** (aka don't repeat the Sepich-Poore et al mistake)

* **MMSeqs2**

* You will need to use MMSeqs2 v13.45111 to match the pipeline

```
conda create -n mmseqsENV mmseqs==13.45111
conda activate mmseqsENV
mmseqs databases UniRef50 <where you want it>/UniRef50 tmp
mmseqs databases UniRef90 <where you want it>/UniRef90 tmp
```

* **CheckM2**

```
conda create -n checkm2ENV checkm2==1.0.1
conda activate checkm2ENV
checkm2 database --download --path <where you want it>
```

* **GTDB-TK**

* Note you can change release214 to release220 (see https://ecogenomics.github.io/GTDBTk/installing/index.html)

```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release214/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
```

* **Bakta**

```
conda create -n baktaENV bakta==1.8.2
conda activate baktaENV
bakta_db download --output <output-path> --type full
```

* **SingleM**

```
conda create -n singlemENV singlem==0.18.3
conda activate singlemENV
singlem data --output-directory <where you want it>
```

# Pipeline

1. With all BAM files (.bam) files as downloaded from the TCGA/CPTAC/where-ever you want in the `input_bams` directory, this will extract all unaligned paired end reads (as marked in the BAM) and generate total read counts for each BAM. The unaligned reads will be in the `UNALIGNED_FASTQ` directory.

```
tcga_cptac_taxonomic_profiler extract --input input_bams --output TCGA_output 
```

2. Download the CHM13 `human-t2t-hla.fa` host genome (from [hostile](https://github.com/bede/hostile)) and combine it with phix174 (a common contaminant spike-in in Illumina sequencing runs) with the following command.

```
tcga_cptac_taxonomic_profiler install-host --database host_genome_db
```

3. Run [trimnami](https://github.com/beardymcjohnface/Trimnami) specifying the directory of FASTQ reads as `--reads`.

Note that you will need to modify the config file so that `--length_required 40` under `qc` and `fastp`. This is because a lot of the TCGA reads are 50 bp.

```
trimnami config
# edit the file so that --length_required 40
trimnami run --reads TCGA_output/UNALIGNED_FASTQ --host host_genome_db/human-t2t-hla_phix174.fa fastp --output TCGA_output/trimnami_output --configfile trimnami.config.yaml
```

4. Run the profilers

You will find the host depleted and trimmed reads in `TCGA_output/trimnami_output/fastp`.

To run the Kraken based profiling

```
tcga_cptac_taxonomic_profiler kraken --input TCGA_output/trimnami_output/fastp  --output TCGA_output --fastqc
```

To run the MMseqs2 based profiling

```
tcga_cptac_taxonomic_profiler mmseqs --input TCGA_output/trimnami_output/fastp  --output TCGA_output --fastqc
```

To run the singleM based profiling

* Note: the singleM developers [recommend](https://wwood.github.io/singlem/tools/pipe) using the raw reads not the trimmed/quality control reads, so that is what we will do

```
tcga_cptac_taxonomic_profiler singlem --input TCGA_output/UNALIGNED_FASTQ  --output TCGA_output --fastqc
```

5. Assembly 

* Done with metaspades

```
tcga_cptac_taxonomic_profiler assembly --input tcga_trimnami_output/results/fastp   --output TCGA_output_all --rerun-triggers mtime 
```

6. Binning & Assessment

* The sample-assemblies were binned using VAMB v4.1.3. 
    * VAMB v4.1.3 needs to be in the local env/PATH (conda is behind source for VAMB). It is in the `setup.py` file

* After binning, CheckM2 and GTDB-TK (release 214) are run to assess the quality of the bins

```
tcga_cptac_taxonomic_profiler binning --input tcga_trimnami_output/results/fastp   --output TCGA_output_all 
```


