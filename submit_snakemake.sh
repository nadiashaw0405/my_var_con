#!/bin/bash
#$ -cwd
#$ -j y
#$ -o logs/snakemake_main.log
#$ -l h_rt=24:00:00,h_data=8G,highp
#$ -N varcon_snakemake

# 1. Load Anaconda and activate environment
. /u/local/Modules/default/init/modules.sh
module load anaconda3
source /u/local/apps/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate snakemake_env

# Re-initialize module
export PATH=$PATH:/u/local/bin:/usr/local/bin
. /u/local/Modules/default/init/modules.sh

# 2. Run Snakemake
snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -l h_data=64G,h_rt=12:00:00 -pe shared {threads} -cwd -j y -o logs/cluster/" \
    --jobs 4 \
    --default-resources \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds \
    --verbose \
    --cores 16
