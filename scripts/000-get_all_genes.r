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

# Save in 1000 gene chunks
egenes = unique(sig.sumstat.df$phenotype_id)
chunk_size <- 1000
n_chunks <- ceiling(length(egenes) / chunk_size)

for (i in 1:n_chunks) {
  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, length(egenes))
  chunk <- egenes[start_idx:end_idx]
  
  # Save to file (one gene per row, no header)
  write.table(chunk, 
              file = paste0(outdir, "/all_egenes_chunk_", i, ".txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}

cat("Saved", n_chunks, "files with", length(egenes), "total genes\n")

# Also save a map of gene to chromosome
gene_chr_map <- sig.sumstat.df %>%
  rowwise() %>% 
  mutate(
    chromosome = unlist(strsplit(variant_id, "\\:"))[1],
    chromosome = gsub("chr", "", chromosome)
  ) %>% 
  ungroup() %>%
  select(phenotype_id, chromosome) %>%
  distinct()


write.table(gene_chr_map,
            file = paste0(outdir, "/gene_chr_map.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t")