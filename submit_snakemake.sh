#!/bin/bash
#$ -cwd
#$ -j y
#$ -o logs/snakemake_main.log
#$ -l h_rt=24:00:00,h_data=8G,highp
#$ -N varcon_snakemake

# 1. Load Anaconda and activate your manager environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake_env

# 2. Run Snakemake
snakemake \
    --jobs 4 \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds \
    --cluster-generic-submit-cmd "qsub -l h_data=16G,h_rt=12:00:00 -pe shared {threads} -cwd -j y -o logs/cluster/" \
    --verbose
