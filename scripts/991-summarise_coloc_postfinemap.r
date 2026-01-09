# module load $uber
library(dplyr)
library(tidyverse)
library(ggplot2)

# Options
csdir = "results/coloc"
varex_f = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/plots/multi_tissue_2025/coloc/coloc_loci/colocs_table-var_explained-0pt75.tsv"
genes_f = "input/genes.txt"
eqtl_dir = "../../../../core_analysis_output/IBDverse_multi-tissue_eQTL_project/2025_06_11-multi_tissue_base_results/"
coloc.dir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed/'
vcf.path = '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/imputed.vcf.gz'
repo.dir = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/IBDVerse-sc-eQTL-code"
data.dir = paste0(repo.dir, "/data/")
annotation_file = paste0(data.dir, "all_IBDverse_annotation_mastersheet.csv")
source(paste0(repo.dir,'/qtl_plot/helper_functions.R'))

################
# Load in the coloc and variance explained results (to get the target coloc lead)
################
csfiles = list.files(csdir, pattern="_coloc_postfinemap.txt", full.names=TRUE)
cs_results = lapply(csfiles, read.delim) 
cs_results = Filter(function(x) nrow(x) > 0, cs_results)

varex = read.delim(varex_f) %>%
    mutate(phenotype_id_variant = paste0(phenotype_id, "_", variant_id))
genes = readLines(genes_f)

################
# Summarise fine mapping success (across all, and per resolution)
################
# Look at overall success rate
length(cs_results) # Number of genes that could be fine mapped in at least one condition.
length(genes) # Total number of genes we tried to fine map

# load in the colocs and get list of genes that are qvalue < 0.05 at each resolution
known.coloc.df <- get_colocs(coloc.dir, known_ibd_only = TRUE) %>%
    mutate(gene_level = paste0(phenotype_id, "_", annotation_type))
known.coloc.df.sig = known.coloc.df %>% 
  filter(PP.H4.abf>0.75) # Filter

gene_levels = known.coloc.df.sig %>% 
    filter(phenotype_id %in% genes) %>%
    group_by(phenotype_id) %>%
    slice_min(Level, with_ties=FALSE) %>%
    distinct(phenotype_id, annotation_type, gene_level) # find the minimum resolution at which a gene colocalises with disease

# Merge the cs results, convert celltype names
cs = bind_rows(cs_results) %>%
    left_join(
        annot.mapping, by = c("celltype" = "label_machine")
    ) %>% 
    mutate(
        tissue = case_when(
            grepl("_blood", celltype) ~ "Blood",
            grepl("_r", celltype) ~ "Rectum",
            grepl("_ti", celltype) ~ "TI",
            TRUE ~ "Cross-site"
        ),
        format_name = ifelse(annotation_type == 1, paste0(label_new, " major population"), label_new),
        format_name = paste0(format_name, " (", tissue, ")")
    ) %>%
    mutate(
        gene_level = paste0(phenotype_id, "_", annotation_type)
    ) 

# Filter for credible sets that contain the colocalising variant
cscoloc = cs %>%
    mutate(phenotype_id_variant = paste0(phenotype_id, "_", variant_id)) %>%
    filter(phenotype_id_variant %in% varex$phenotype_id_variant)

length(unique(cscoloc$phenotype_id)) # The number of colocalising eQTLs that could be fine mapped

# Have a look at the number of genes that could be fine mapped at each resolution
cscoloc %>% 
    distinct(phenotype_id, annotation_type) %>%
    group_by(annotation_type) %>%
    summarise(n_finemapped_genes = n())

cs_minres = cscoloc %>%
    filter(gene_level %in% gene_levels$gene_level) # Filter the credible sets for the minimum resolution at which a gene colocalises with disease

fm_gene_levels = cs_minres %>%
    distinct(phenotype_id, annotation_type, gene_level)

# Look at success per resolution
levels = names(table(gene_levels$annotation_type))
fm_sum = data.frame(
    annotation_type = levels,
    n_sig_genes = as.vector(table(gene_levels$annotation_type)),
    n_finemapped_genes = as.vector(table(fm_gene_levels$annotation_type))
) %>% 
    mutate(
        rate = n_finemapped_genes / n_sig_genes
    )

################
# Identify genes with colocalising variants with different directions of effect
################
coloc_leads = varex %>% 
    pull(phenotype_id_variant) %>% unique()

# For each gene, first check the direction of effect of the hash that contains the exact colocalising variant 
cs_doe_sum = cscoloc %>% group_by(phenotype_id) %>% 
    mutate(phenotype_id_variant = paste0(phenotype_id, "_", variant_id)) %>%
    filter(phenotype_id_variant %in% coloc_leads) %>% 
    group_by(phenotype_id, strict_collapsed_hash) %>%
    mutate(direction_of_effect = ifelse(slope > 0, "positive", "negative")) %>%
    summarise(
        npositive = sum(direction_of_effect == "positive"),
        nnegative = sum(direction_of_effect == "negative"),
        positive_celltypes = paste(unique(celltype[direction_of_effect == "positive"]), collapse=","),
        negative_celltypes = paste(unique(celltype[direction_of_effect == "negative"]), collapse=",")
    ) 

opposite_coloc_variant = cs_doe_sum %>%
    filter(npositive > 0 & nnegative > 0) 

opposite_coloc_variant %>% pull(phenotype_id) %>% unique() %>% length() # 5 - This is the number of disease effector genes for which we do see opposite direction of effect for the colocalising variant AND the signals are fine mapped, AND they colocalise.

################
# Deep dive into specific gene
################
gene="ENSG00000115604"
gene_name="IL18R1"
coloc_variant = coloc_leads[grep(gene, coloc_leads)] %>% strsplit(., "\\_") %>% unlist() %>% tail(1)
disease = varex %>% filter(phenotype_id == !!gene, variant_id == !!coloc_variant) %>% pull(gwas_trait) %>% unique()
coloc_disease = known.coloc.df.sig %>% mutate(qtl_lead = gsub("\\_", "\\:", qtl_lead)) %>% filter(phenotype_id == !!gene, qtl_lead == !!coloc_variant, gwas_trait == !!disease) %>% slice_max(PP.H4.abf) %>% pull(PP.H4.abf)
intcs  = cscoloc %>% 
    filter(
        strict_collapsed_hash == opposite_coloc_variant %>% filter(phenotype_id == !!gene) %>% pull(strict_collapsed_hash),
        phenotype_id == !!gene,
        variant_id == !!coloc_variant
    )


# ---- 1. Check the coloc between these two conditions
topup = intcs %>% 
    filter(slope > 0) %>% 
    slice_max(PIP) %>% 
    pull(celltype)

topdown = intcs %>% 
    filter(slope < 0) %>% 
    slice_max(PIP) %>% 
    pull(celltype)

colocres = read.delim(paste0(csdir, "/", gene, "_pairwise_coloc_results.txt")) %>%
    filter(
        (cond1 == !!topup & cond2 == !!topdown) | (cond1 == !!topdown & cond2 == !!topup)
    )

print(paste0("..The coloc between ", topup, " and ", topdown, " for variant ", coloc_variant, " is PP.H4=", round(colocres$PP.H4.abf, 3)))
print(paste0("..", gene_name, " colocalises with ", disease, " with PP.H4=", round(coloc_disease,3)))

# ---- 2. Plot the eQTL from this variant in the interesting conditions

# Read in the expression for just this gene, and merge with genotypes for just this variant
exprlist = lapply(c(topup, topdown), function(ct) {
    # Expression
    exprpath = paste0(eqtl_dir, "/dMean__", ct, "_all/OPTIM_pcs/base_output/base/Expression_Data.sorted.bed")
    exprtemp = paste0("temp/", gene, "-", ct, "_expression.tsv")
    system(sprintf("head -1 %s > %s", exprpath, exprtemp))
    system(sprintf("grep '%s' %s >> %s", gene, exprpath, exprtemp))
    subexpr = read.delim(exprtemp) %>% 
        select(-c("X.chr", "start", "end")) %>% 
        tidyr::pivot_longer(cols = -gene_id, names_to = "genotype_id", values_to = "expression") %>%
        mutate(
            genotype_id = gsub("X3", "3", genotype_id),
            condition = ct
        )

    # Genotypes
    spec_var_pos = paste(unlist(strsplit(as.character(coloc_variant), "\\:"))[c(T,T,F,F)], collapse=":")
    specific_vcf_path = paste0("temp/", coloc_variant, ".vcf.gz")
    if(!file.exists(specific_vcf_path)){
        system(sprintf('bcftools view -r %s %s -Oz -o %s', spec_var_pos, vcf.path, specific_vcf_path))
    }

    # Load vcf
    vcf.df <- read.table(specific_vcf_path,skip=50, header=TRUE, comment.char = "")
    
    variant.df <- filter(vcf.df, ID == coloc_variant) %>% 
        pivot_longer(!names(vcf.df)[1:9] , names_to = 'genotype_id') %>% 
        mutate(genotype_id = gsub('^X', '', genotype_id)) %>% 
        separate(value, into = str_split('GT:DS:HDS:GP', ':')[[1]],sep = ':') %>% 
        mutate(Genotype = str_count(GT, "1")) %>% 
        select(genotype_id,Genotype) %>% 
        drop_na() 
    
    # merge
    merged.df <- left_join(subexpr, variant.df, by = "genotype_id")
    return(merged.df)

})
expr = bind_rows(exprlist) %>% 
    filter(!is.na(Genotype))

# Annotate the condition with proper name and category
expr = expr %>%
    left_join(annot.mapping, by = c("condition" = "label_machine")) %>% 
    mutate(
        tissue = case_when(
            grepl("_blood", condition) ~ "Blood",
            grepl("_r", condition) ~ "Rectum",
            grepl("_ti", condition) ~ "TI",
            TRUE ~ "Cross-site"
        ),
        format_name = ifelse(annotation_type == 1, paste0(label_new, " major population"), label_new),
        format_name = paste0(label_new, " (", tissue, ")")
    )

# Plot a boxplot, like the variance explained
# Prepare annotations with slope and p-values
annotations = intcs %>%
    filter(celltype %in% expr$condition) %>%
    mutate(
        annotation_text = sprintf("β = %.3f, P = %.2e", slope, pval_nominal)
    ) %>%
    select(format_name, annotation_text)

# Filter expr to only include format_names present in the data
expr = expr %>%
    filter(format_name %in% annotations$format_name) %>% 
    mutate(format_name = as.character(format_name))

outf=paste0("results/plots/", gene, "-", gene_name, "_opposite_coloc_DoE.png")
ggplot(expr, aes(x = factor(Genotype), y = expression, fill = category_new)) +
    geom_boxplot() + 
    facet_grid(. ~ format_name, scales = "free_x", space = "free_x") + 
    geom_text(data = annotations, aes(x = 2, y = Inf, label = annotation_text), 
              inherit.aes = FALSE, vjust = 1.5, size = 5, fontface = "plain") +
    theme_classic() + 
    labs(
        title = bquote("Strongest colocalising, opposing effects on"~italic(.(gene_name))),
        x = paste0("Genotype at colocalising variant (", coloc_variant, ")"),
        y = bquote("Rank-normalised pseudobulked"~italic(.(gene_name))~"expression")
    ) + 
    scale_fill_manual(values = umap.category.palette) + 
    theme(
        legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 15)
    )

ggsave(outf, width = 7, height = 6)

# ---- 3. Plot the credible sets from these two conditions
outf=paste0("results/plots/", gene, "-", gene_name, "_opposite_coloc_Manhattan.png")
coloc_hash = intcs %>% pull(strict_collapsed_hash) %>% unique()
plt_cs = cs %>% 
    filter(
        strict_collapsed_hash %in% !!coloc_hash,
        celltype %in% c(topup, topdown)
    ) %>% 
    mutate(
        coloc_qtl = variant_id == coloc_variant,
        position = unlist(strsplit(as.character(variant_id), "\\:"))[c(F,T,F,F)] %>% as.numeric()
    ) 
    

ggplot(plt_cs, aes(x=position, y=PIP, color=coloc_qtl)) +
    geom_point(size=3) +
    facet_grid(format_name~.) + 
    theme_classic() +
    scale_color_manual(values=c("grey", "darkorange2")) +
    labs(
        title = bquote("Regional association of colocalising eQTLs for"~italic(.(gene_name))~""),
        x = "Genomic position",
        y = "PIP"
    ) +
    theme(
        legend.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 15)
    ) + 
    theme(
        legend.position = "none"
    )

ggsave(outf, width = 7, height = 6)


################
# Deep dive into FUBP1
################
gene="ENSG00000162613"
gene_name="FUBP1"
coloc_variant = coloc_leads[grep(gene, coloc_leads)] %>% strsplit(., "\\_") %>% unlist() %>% tail(1)
disease = varex %>% filter(phenotype_id == !!gene, variant_id == !!coloc_variant) %>% pull(gwas_trait) %>% unique()
coloc_disease = known.coloc.df.sig %>% mutate(qtl_lead = gsub("\\_", "\\:", qtl_lead)) %>% filter(phenotype_id == !!gene, qtl_lead == !!coloc_variant, gwas_trait == !!disease) %>% slice_max(PP.H4.abf) %>% pull(PP.H4.abf)

# How many conditions can FUBP1 eQTLs be found in?
cs %>% filter(phenotype_id == !!gene) %>% distinct(celltype) %>% nrow() # 11

# How many fine mapped eQTLs?
cs %>% 
    filter(
        phenotype_id == !!gene
    ) %>% 
    distinct(strict_collapsed_hash) %>% nrow() # 16

# How many fine mapped eQTLs that contain the coloc variant?
cs %>%
    filter(
        phenotype_id == !!gene,
        variant_id == !!coloc_variant
    ) %>%
    distinct(strict_collapsed_hash) %>% nrow() # 2. However, this is because of one set that hasn't merged properly (see below. So overall, think there are 15 eQTLs, 1 with colocalising eQTL)

# What is the range of PIPs for the colocalising variant across conditions?
cs %>%
    filter(
        phenotype_id == !!gene
    ) %>%
    group_by(celltype, credible_set) %>%
    mutate(
        rank = rank(-PIP)
    ) %>%
    filter(
        variant_id == !!coloc_variant
    ) %>%
    select(celltype, PIP, rank)  

# Look at the finemapped sets that contain the disease-associated variant
intcs = cscoloc %>% 
    filter(
        strict_collapsed_hash %in% c(cs_doe_sum %>% filter(phenotype_id == !!gene) %>% pull(strict_collapsed_hash)),
        phenotype_id == !!gene,
        variant_id == !!coloc_variant
    ) 

# The merging here is a bit complicated. the colocalising variant for FUBP1 seems to be a part of two collapsed hash's. Check the coloc between these
colocres = read.delim(paste0(csdir, "/", gene, "_pairwise_coloc_results.txt")) %>%
    filter(
        hash1 %in% intcs$hash & hash2 %in% intcs$hash | hash2 %in% intcs$hash & hash1 %in% intcs$hash
    )

# It looks like something has gone wrong with the merging here, combine these hashes manually
intcs$collapsed_hash = paste0(unique(intcs$collapsed_hash), collapse="|")
intcs$strict_collapsed_hash = paste0(unique(intcs$collapsed_hash), collapse="|")