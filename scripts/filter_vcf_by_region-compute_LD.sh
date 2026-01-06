#!/bin/bash

# Filter VCF file for specific genomic regions
# Usage: ./filter_vcf_by_region.sh <input_vcf> <region_file> <output_vcf>
# Region file format: tab-separated, no header (chromosome start end)

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_vcf> <region_file> <output_vcf>"
    echo ""
    echo "Region file format: tab-separated, no header (chromosome start end)"
    exit 1
fi

INPUT_VCF="$1"
REGION_FILE="$2"
OUTPUT_VCF_PREFIX="$3"
STATS_VCF_FILE="${OUTPUT_VCF_PREFIX}-vcfstats_prefiltering.txt"
STATS_VCF_OUTPUT_FILE="${OUTPUT_VCF_PREFIX}-vcfstats_postfiltering.txt"
LDMATRIX_FNAME="$4"

# Check if files exist
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF file '$INPUT_VCF' not found"
    exit 1
fi

if [ ! -f "$REGION_FILE" ]; then
    echo "Error: Region file '$REGION_FILE' not found"
    exit 1
fi


# Filter VCF file using region strings (convert to comma-separated format)
# bcftools -r expects comma-separated regions in chr:start-end format
REGIONS=$(awk '{printf "%s:%s-%s,", $1, $2, $3}' "$REGION_FILE" | sed 's/,$//')

echo "Filtering regions:"
echo "$REGIONS"

bcftools view -r "$REGIONS" "$INPUT_VCF" -o "${OUTPUT_VCF_PREFIX}_filtered.vcf"

if [ $? -eq 0 ]; then
    echo "Successfully filtered VCF!"
    echo "..Output: ${OUTPUT_VCF_PREFIX}_filtered.vcf"
else
    echo "Error: VCF filtering failed"
    exit 1
fi

# Then filter the vcf for low frequency variants
echo ".. Filtering for low frequency variants (MAF < 0.005)"
bcftools view -q 0.005:minor "${OUTPUT_VCF_PREFIX}_filtered.vcf" -o "${OUTPUT_VCF_PREFIX}_filtered-low-mafs.vcf.gz"

bcftools stats ${OUTPUT_VCF_PREFIX}_filtered.vcf > "$STATS_VCF_FILE"
bcftools stats "${OUTPUT_VCF_PREFIX}_filtered-low-mafs.vcf.gz" > "$STATS_VCF_OUTPUT_FILE"

echo "Successfully filtered VCF!"
echo "..Output: $VCF_OUTPUT_FILE"
echo "..Stats files: $STATS_VCF_FILE, $STATS_VCF_OUTPUT_FILE"

# Convert the vcf to bgen
echo "..Converting to bgen format"
module load HGI/softpack/users/cc53/finemapping/1
qctool -g "${OUTPUT_VCF_PREFIX}_filtered-low-mafs.vcf.gz" -og "${OUTPUT_VCF_PREFIX}_genotypes.bgen" -os "${OUTPUT_VCF_PREFIX}_genotypes.sample"
bgenix -index -g "${OUTPUT_VCF_PREFIX}_genotypes.bgen" -clobber
echo "..BGEN files created: ${OUTPUT_VCF_PREFIX}_genotypes.bgen, ${OUTPUT_VCF_PREFIX}_genotypes.sample"
echo "..VCF filtering and conversion to BGEN completed!"

# get a plink2 .bim file to make the .z input file for LDstore2
echo "..Creating PLINK2 .bim file"
plink2 --bgen ${OUTPUT_VCF_PREFIX}_genotypes.bgen ref-unknown --make-just-bim --out ${OUTPUT_VCF_PREFIX}

#make the .z input file
echo "..Making .z input file"
echo "rsid chromosome position allele1 allele2" > ${OUTPUT_VCF_PREFIX}.z
awk '{print $2,"chr"$1,$4,$5,$6}' ${OUTPUT_VCF_PREFIX}.bim >> ${OUTPUT_VCF_PREFIX}.z

#make the master file
echo "..Creating master file for LDstore2"
echo "z;bgen;bgi;bcor;ld;n_samples" > ${OUTPUT_VCF_PREFIX}_masterfile.txt
echo "$OUTPUT_VCF_PREFIX.z;${OUTPUT_VCF_PREFIX}_genotypes.bgen;${OUTPUT_VCF_PREFIX}_genotypes.bgen.bgi;$LDMATRIX_FNAME.bcor;$LDMATRIX_FNAME.ld;456" >> ${OUTPUT_VCF_PREFIX}_masterfile.txt

#execute LDstore2
echo "..Executing LDstore2"
ldstore2 --in-files ${OUTPUT_VCF_PREFIX}_masterfile.txt --read-only-bgen --n-threads 1 --write-bcor
ldstore2 --in-files ${OUTPUT_VCF_PREFIX}_masterfile.txt --bcor-to-text
echo "..LDstore2 computation completed!"


