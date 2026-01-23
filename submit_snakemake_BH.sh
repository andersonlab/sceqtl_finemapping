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
config_var=configs/config_sigpairs.yaml
worfklow_prefix="FM_"
group="team152"
workdir=${PWD}

# Load snakemake and singulatiry
module load HGI/common/snakemake/7
module load ISG/singularity/3.11.4
which singularity

# Detect the nature of the analysis to be done. 
# If config has the option 'gene_cond_file', then we are testing gene x condition pairs explicitly, not all vs all.
if grep -q "gene_cond_file" ${config_var}; then
    workflow_snakefile="workflow/Snakefile_sigpairs"
    gene_cond_file=$(grep "gene_cond_file" ${config_var} | awk '{print $2}' | tr -d '"')
    npairs=$(zcat ${gene_cond_file} | wc -l)
    echo "-- Running pipeline on ${npairs} gene x condition pairs from config file. --"
else
    workflow_snakefile="workflow/Snakefile"
    echo "-- Running snakemake with all genes and all conditions. --"
fi

# Make a log dir
mkdir -p sm_logs

# Run snakemake
snakemake -j 5000 \
    --latency-wait 90 \
    --use-envmodules \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --keep-going \
    --directory ${workdir} \
    --cluster-config ${config_var} \
    --cluster-config configs/cluster_config.yaml \
    --use-singularity \
    --singularity-args "-B /lustre,/software" \
    --restart-times 3 \
    --snakefile ${workflow_snakefile}

# bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q basement -J "snakemake_master_FINEMAP" < submit_snakemake_BH.sh 
