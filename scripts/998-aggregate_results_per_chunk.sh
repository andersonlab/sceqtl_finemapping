#!/bin/bash
# Aggregate all results per chunk into single files (also removing empty files)

RESDIR="$1" # RESDIR=results_Naturepaper_copy
CHROM_MAP="$2" # CHROM_MAP=input/gene_chr_map.txt 

# For each gene, aggregate all the finemapping variants into a single file
echo "..Aggregating results per gene"
for gene in $RESDIR/finemapping/*; do
    gene_name=$(basename ${gene})
    echo "....Processing gene: ${gene_name}"
    # Aggregate all credible sets into a single file
    output_file=${RESDIR}/finemapping/${gene_name}/${gene_name}_all_crediblesets.txt
    rm -f ${output_file}  # Remove existing file if any
    first_file=true 
    for cs_file in ${RESDIR}/finemapping/${gene_name}/*crediblesets.txt; do
        if [ -f "${cs_file}" ]; then
            if [ "$first_file" = true ]; then   
                # Get header from first file and add 'condition' column
                head -n 1 "${cs_file}" | awk '{print $0"\tcondition"}' > ${output_file}
                first_file=false
            fi
            # Extract condition name from filename (everything before first '.')
            condition=$(basename "${cs_file}" | cut -d'.' -f1)
            tail -n +2 "${cs_file}" | awk -v cond="$condition" '{print $0"\t"cond}' >> ${output_file}  # Skip header and add condition
        fi
    done
done

# Remove all the individual credible set files after aggregation
echo "..Removing per-cell-type credible set files"
find ${RESDIR}/finemapping/*/ -maxdepth 1 -name "*.crediblesets.txt" -exec rm {} \;




