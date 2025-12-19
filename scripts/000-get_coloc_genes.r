library(tidyverse)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

variance_explained_f = args[1] # variance_explained_f = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv"

# Load
varex = read.delim(variance_explained_f)

# Find genes with opposing directions of effect
test_genes  = varex %>% 
    group_by(phenotype_id) %>% 
    filter(!(all(beta > 0) | all(beta < 0))) %>%
    pull(phenotype_id) %>% 
    unique() 

write.table(test_genes, "input/genes.txt", sep = "\n", quote=F, row.names=F, col.names=F)
