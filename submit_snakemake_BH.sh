#!/bin/bash
#BSUB -o sm_logs/snakemake_master-%J-output.log
#BSUB -e sm_logs/snakemake_master-%J-error.log 
#BSUB -q oversubscribed
#BSUB -G team152
#BSUB -n 1
#BSUB -M 10000
#BSUB -a "memlimit=True"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -J 1

# Define some params
config_var=config.yaml
worfklow_prefix="FM_"
group="team152"
workdir=${PWD}

# Load snakemake and singulatiry
module load HGI/common/snakemake/7
module load ISG/singularity/3.11.4
which singularity


# Make a log dir
mkdir -p sm_logs

# Run snakemake
snakemake -j 20000 \
    --latency-wait 90 \
    --use-envmodules \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --keep-going \
    --directory ${workdir} \
    --cluster-config ${config_var} \
    --cluster-config cluster_config.yaml \
    --use-singularity \
    --singularity-args "-B /lustre,/software" \
    --keep-going \
    --restart-times 0 

# bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_FINEMAP" < submit_snakemake_BH.sh 
