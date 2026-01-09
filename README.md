### Fine mapping sc-eQTLs and colocalisation of credible sets.
#### Bradley - 09/01/2026
Adjustment of Celeste E. Cohens's pipeline to be ran with snakemake and on specific genes only across all cell-types.
Requires eQTLs to already be calculated.

### Installation
If working on Sanger farm, no instals needed. Otherwise, need to install the finemapping singularity container from docker and adjust paths in the snakefile accordingly. 
```
export SINGULARITY_CACHEDIR=$PWD/.singularity_cache
mkdir -p "$SINGULARITY_CACHEDIR"
singularity pull docker://bh18/sceqtl_finemapping
```

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

