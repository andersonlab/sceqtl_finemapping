### Bradley
Adjustment of Celeste's pipeline - fine mapping specific gene in all cell-types (regardless of significance). 
1. Define celltypes
```
mkdir -p {logs,results,input}
eqtl_dir="../../../../core_analysis_output/IBDverse_multi-tissue_eQTL_project/2025_06_11-multi_tissue_base_results"
for f in ${eqtl_dir}/*; do
  f=$(basename "${f#dMean__}")
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




1. Subset vcf to region around specific gene
```
module load $uber
mkdir -p input/{vcf,LDmatrix}
mkdir -p results/finemapping
mkdir -p logs
Rscript scripts/get_gene_window.r $gene $gene_ens
```

2. Filter vcf for the window, remove low MAF variants, convert to bgenix, compute LD
```
module load bcftools-1.19
filt_vcf="input/vcf/${gene}_"
LD_file_pref="input/LDmatrix/${gene}"
region_file="input/vcf/${gene}_gene_window.txt"
bash scripts/filter_vcf_by_region-compute_LD.sh "$VCF_FILE" "$region_file" "$filt_vcf" "$LD_file_pref"
```

3. Finemap
```
MEM=5000
fmpath="results/finemapping"
mkdir -p ${fmpath}/$gene
module load HGI/softpack/users/cc53/finemapping_R/2
for celltype in "${format_celltypes[@]}"; do
    echo ".. Submitting finemapping of $gene in $celltype"
    SAMPLESIZE=$(awk '{print NF-1; exit}' $eqtl_dir/dMean__${celltype}_all/OPTIM_pcs/base_output/base/Covariates.tsv)
    chr=$(awk -F'\t' 'NR==1{ sub(/^chr/,"",$1); print $1 }' input/vcf/${gene}_gene_window.txt)
    bsub -J "finemap-${gene}-${celltype}" -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" -G team152 \
    -e logs/finemap-${gene}-${celltype}-stderr \
    -o logs/finemap-${gene}-${celltype}-stdout \
        "Rscript scripts/runsusie.r \
                    $gene_ens \
                    $eqtl_dir/dMean__${celltype}_all/OPTIM_pcs/base_output/base/cis_nominal1.cis_qtl_pairs.$chr.tsv \
                    $SAMPLESIZE \
                    ${LD_file_pref}.ld \
                    input/vcf/${gene}_.z \
                    ${fmpath}/${gene}/${celltype}."
done
```

