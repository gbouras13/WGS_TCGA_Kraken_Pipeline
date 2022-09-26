#!/bin/bash -l

#SBATCH --job-name=coad_tcga_snkg
#SBATCH --mail-user=george.bouras@adelaide.edu.au
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --err="coad_tcga_snk.err"
#SBATCH --output="coad_tcga_snk.out"

# Resources allocation request parameters
#SBATCH -p batch
#SBATCH -N 1              	                                # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 1              	                                # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=2-23:59:50                                         # time allocation, which has the format (D-HH:MM), here set to 1 hou                                           # generic resource required (here requires 1 GPUs)
#SBATCH --mem=1GB                                              # specify memory required per node


SNK_DIR="/hpcfs/users/a1667917/Kevin/WGS_TCGA_Kraken_Pipeline"
PROF_DIR="/hpcfs/users/a1667917/snakemake_slurm_profile"

cd $SNK_DIR

module load Anaconda3/2020.07
conda activate snakemake_clean_env

snakemake  -c 1 -s /hpcfs/users/a1667917/Kevin/WGS_TCGA_Kraken_Pipeline/extract_unaligned_fastq.smk \
--config Bams=/hpcfs/users/a1667917/Kevin/COAD_Bams/All_Bams/ Output=/hpcfs/users/a1667917/Kevin/COAD_WGS_Output_Final \
--profile $PROF_DIR/wgs_tcga \
--keep-going

# only need the output directory after extract_unaligned_fastq.smk has been run

snakemake  -c 1 -s /hpcfs/users/a1667917/Kevin/WGS_TCGA_Kraken_Pipeline/wgs_runner.smk \
--config Output=/hpcfs/users/a1667917/Kevin/COAD_WGS_Output_Final \
--profile $PROF_DIR/wgs_tcga \
--keep-going

conda deactivate
