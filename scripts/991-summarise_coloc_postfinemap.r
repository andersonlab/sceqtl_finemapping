library(dplyr)

# Options
csdir = "results/coloc"
varex_f = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv"
genes_f = "input/genes.txt"

################
# Load in the coloc and variance explained results (to get the target coloc lead)
################
csfiles = list.files(csdir, pattern="_coloc_postfinemap.txt", full.names=TRUE)
cs_results = lapply(csfiles, read.delim)

varex = read.delim(varex_f)
genes = readLines(genes_f)

################
# Summarise fine mapping success (across all, and per resolution)
################
length(cs_results) # Number of genes that could be fine mapped in at least one condition.

# Look at success per resolution

################
# Identify genes with colocalising variants with different directions of effect
################
cs = do.call(rbind, cs_results)
# For each gene, get the hash that contains the colocalising variant



# Then see whether the variants from that credible set hash (or just the lead) have the same or different directions of effect across conditions.