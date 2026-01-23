### Fine mapping sc-eQTLs and colocalisation of credible sets.
#### Bradley - 09/01/2026
Adjustment of Celeste E. Cohens's pipeline to be ran with snakemake and on specific genes only across all cell-types.
Requires eQTLs to already be calculated. \
\
*** WARNING: Depending on the size of the cohort, per-gene temporary LD files (`.ld`) are often very large (~0.5Gb for a cohort of ~400 individuals). While as many temporary files are removed during running of the pipeline as possible, this means that if you have many genes being tested the working directory can require a substantial amount of memory. Please ensure there is enough quota to complete, else the pipeline will fail. For 20k genes, 250 cell-types and 400 individuals, this requires ~ 7TB free memory. *** 

### Installation
If working on Sanger farm, no installs needed. Otherwise, need to install the finemapping singularity container from docker and adjust paths in the snakefile accordingly. 
```
export SINGULARITY_CACHEDIR=$PWD/.singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"
singularity pull docker pull bh18/sceqtl_finemapping:5
```

### If running on specific sets of genes (e.g hundreds of disease effector genes)
1. Define celltypes
```
mkdir -p {logs,results,input}
eqtl_dir="../../../../core_analysis_output/IBDverse_multi-tissue_eQTL_project/2025_06_11-multi_tissue_base_results"
for f in ${eqtl_dir}/*; do
  f=$(basename $f)
  f=${f#dMean__}
  f="${f%_all}"
  echo $f >> input/celltypes.txt
done
```

2. Get the genes we want to test
```
varex_f="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv"
Rscript scripts/000-get_coloc_genes.r $varex_f
```

3. Submit pipeline
```
bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_FINEMAP" < submit_snakemake_BH.sh 
```

### If running on all genes - need to chunk
1. Define cell-types as before. Then define genes using all eGenes in 1000 gene chunks. These files are `input/all_egenes_chunk_<i>.txt`
```
Rscript scripts/000-get_all_genes.r
```

2. Submit pipeline. NOTE: Need to adjust the gene_file flag in `config.yaml` to run for each chunk.
```
bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_FINEMAP" < submit_snakemake_BH.sh 
```

3. Tidy up the output after each chunk.
```
chunk=1
mkdir -p results_all_egenes
mkdir -p results_all_egenes/{coloc,finemapping}
rm -r logs/* # Remove log files
rm -r sm_logs/*
rm -r input/LDmatrix # Remove temporary LD files
rm -r input/vcf # Remove temporary genotyping files
mv results results_${chunk} # Rename
bash scripts/998-aggregate_results_per_chunk.sh results_${chunk} input/gene_chr_map.txt # Aggregate the per-cell-type finemapping files from this chunk and aggregate
rsync -av results_${chunk}/coloc/ results_all_egenes/coloc/ # move results
rsync -av results_${chunk}/finemapping/ results_all_egenes/finemapping/
rm -r results_${chunk} # Remove this chunks data
```

4. Combine results after ALL chunks
```
```

### If running on all genes - significant only
1. Get the significant gene x cell-type pairs
```
mkdir -p input
Rscript scripts/000-get_eGene_cond_pairs.r
```

2. Then run the pipeline (will specify which Snakefile to use, need to just point to the output of the above within `configs/config_sigpairs.yaml`)
```
bsub -M 10000 -a "memlimit=True" -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o sm_logs/snakemake_master-%J-output.log -e sm_logs/snakemake_master-%J-error.log -q oversubscribed -J "snakemake_master_FINEMAP" < submit_snakemake_BH.sh 
```






