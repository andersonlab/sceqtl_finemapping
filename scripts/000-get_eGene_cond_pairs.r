# Bradley Jan 2026
# Set up
library(dplyr)
repo.dir <- '../IBDVerse-sc-eQTL-code/'
sumstats.all.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
data.dir <- paste0(repo.dir,'/data/')
source(paste0(repo.dir,'qtl_plot/helper_functions.R'))
outdir = "input"

# Load eQTLs
sumstat.df = read_eqtls(sumstats.all.basedir, quick=T)
sig.sumstat.df <- sumstat.df %>% 
  filter(qval < 0.05) 

# Save a list of sig gene x condition pairs
gene_cond_pairs <- sig.sumstat.df %>%
    select(phenotype_id, annotation) %>%
    distinct() %>% 
    write.table(
        file = gzfile(paste0(outdir, "/sig_gene_condition_pairs.txt.gz")),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )

# Also a 10 commonly detected eGenes for testing
common_egenes <- sig.sumstat.df %>%
    count(phenotype_id) %>%
    top_n(10, n) %>%
    pull(phenotype_id) 

uncommon <- sig.sumstat.df %>%
    count(phenotype_id) %>%
    arrange(n) %>%
    head(10) %>%
    pull(phenotype_id) 

sig.sumstat.df %>%
    select(phenotype_id, annotation) %>%
    distinct() %>% 
    filter(phenotype_id %in% c(common_egenes, uncommon)) %>% 
    write.table(
        file = gzfile(paste0(outdir, "/common_uncommon_egenes.txt.gz")),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )