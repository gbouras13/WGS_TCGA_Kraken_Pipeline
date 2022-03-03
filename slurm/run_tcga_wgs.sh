#!/bin/bash -l

#SBATCH --job-name=wgs_tcga_snkg
#SBATCH --mail-user=george.bouras@adelaide.edu.au
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --err="wgs_tcga_snk.err"
#SBATCH --output="wgs_tcga_snk.out"

# Resources allocation request parameters
#SBATCH -p batch
#SBATCH -N 1              	                                # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 1              	                                # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=12:00:00                                         # time allocation, which has the format (D-HH:MM), here set to 1 hou                                           # generic resource required (here requires 1 GPUs)
#SBATCH --mem=2GB                                              # specify memory required per node


SNK_DIR="/hpcfs/users/a1667917/Kevin/WGS_TCGA_Kraken_Pipeline"
PROF_DIR="/hpcfs/users/a1667917/snakemake_slurm_profile"

cd $SNK_DIR

module load Anaconda3/2020.07
conda activate snakemake_clean_env

snakemake  -s /hpcfs/users/a1667917/Kevin/WGS_TCGA_Kraken_Pipeline/wgs_runner.smk  --config Reads=/hpcfs/users/a1667917/Kevin/TCGA_WGS_Total_Bams/ Output=/hpcfs/users/a1667917/Kevin/Test_TCGA_WGS_Output  --profile wgs_tcga

conda deactivate
