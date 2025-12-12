### Bradley
Adjustment of Celeste's pipeline - fine mapping specific gene in all cell-types (regardless of significance). Define options
```
gene=FUBP1
gene_ens=ENSG00000162613
VCF_FILE="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/imputed.vcf.gz"
eqtl_dir="../../../../core_analysis_output/IBDverse_multi-tissue_eQTL_project/2025_06_11-multi_tissue_base_results"
```

Get celltypes
```
mapfile -t celltypes < <(ls -1A "$eqtl_dir")
for f in "${celltypes[@]}"; do
  new="${f//dMean__/}"   # remove all occurrences of dMean__
  new="${new//_all/}"    # remove all occurrences of _all
  format_celltypes+=("$new")
done
```

1. Subset vcf to region around specific gene
```
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

