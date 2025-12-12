suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

# Define options
args = commandArgs(trailingOnly=TRUE)
gene = args[1] #gene = "FUBP1"
gene_ens = args[2] # gene_ens = "ENSG00000162613"
window = 1e6
gene_pos_file = "../IBDVerse-sc-eQTL-code/data/gene_counts_Ensembl_105_phenotype_metadata.annotation_file.txt"

# Get coords
read.delim(gene_pos_file) %>% 
    filter(feature_id == gene_ens) %>%
    mutate(window_start = start - window,
           window_end = end + window,
           chromosome = paste0("chr", chromosome)) %>%
    select(chromosome, window_start, window_end) %>% 
    write.table(., file = paste0("input/vcf/", gene, "_gene_window.txt"),
                col.names = F, row.names = FALSE, sep = "\t", quote = FALSE)
    