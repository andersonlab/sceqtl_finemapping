# module load $uber
library(dplyr)
library(tidyr)
library(ggplot2)

# Define gene
gene_of_interest <- "FUBP1"
conditions_file <- "../IBDVerse-sc-eQTL-code/temp/FUBP1_finemap_conditions.txt"
coloc_qtl_variant = "chr1:77984833:C:A"
thresh = "pval_beta,0.05"
nomsumstat_file = paste0("../IBDVerse-sc-eQTL-code/temp/", gene_of_interest, "-", coloc_qtl_variant, "-nominal.txt.gz")

threshcol = unlist(strsplit(thresh, ","))[1]
threshval = as.numeric(unlist(strsplit(thresh, ","))[2])

######## 
# Load in the conditions to test, and the finemapping results
######## 
conditions = read.delim(conditions_file) %>% 
    filter(!!sym(threshcol) < threshval) %>% 
    pull(name)

conditions = gsub("dMean__", "", conditions)
conditions = gsub("_all", "", conditions)
fm = do.call(rbind, lapply(conditions, function(cond) {
    finemap_file = paste0("results/finemapping/", gene_of_interest, "/", cond, ".crediblesets.txt")
    if(file.exists(finemap_file)){
        read.table(finemap_file, header=TRUE, sep=" ", stringsAsFactors=FALSE) %>%
            mutate(condition = cond)
    }
}))

# Count the number of variants in a credible set from each celltype
nvars_per_set = fm %>% 
    group_by(condition, credible_set) %>%
    summarise(nvars = n()) %>%
    ungroup()

fm = fm %>%
    left_join(nvars_per_set, by=c("condition", "credible_set"))

########
# High level summaries
########
nconditions_test = length(conditions) # 67
nconditions_finemapped = length(unique(fm$condition)) # 21 (31%)
nvariants = length(unique(fm$variant_id)) # 36
nvariants_pippt5 = fm %>% filter(PIP > 0.5) %>% pull(variant_id) %>% unique() %>% length() # 2
cs = fm %>%
    group_by(condition, credible_set) %>%
    arrange(variant_id) %>%
    summarise(cs = paste(variant_id, collapse = ",")) %>%
    ungroup() 

ncrediblesets = length(unique(cs$cs)) # 7

topvariant = fm %>%
    group_by(variant_id) %>%
    summarise(
        ncs = n(),
        max_pp = max(PIP),
        median_pp = median(PIP),
        nsolocs = sum(nvars == 1)
    ) %>%
    ungroup() %>% 
    arrange(-nsolocs)

## A tibble: 36 × 5
#   variant_id          ncs max_pp median_pp nsolocs
#   <chr>             <int>  <dbl>     <dbl>   <int>
# 1 chr1:77984833:C:A    19 1.00      0.962       12
# 2 chr1:76981862:G:C     1 0.0644    0.0644       0
# 3 chr1:76982243:C:T     1 0.0647    0.0647       0
# 4 chr1:76983534:A:G     1 0.0648    0.0648       0

# Write out summaries
topvariant %>% slice_max(ncs) %>% pull(variant_id) # Most common: chr1:77984833:C:A (19/21 credible sets)
topvariant %>% slice_max(nsolocs) %>% pull(variant_id) # Most solitary: chr1:77984833:C:A (12/19 credible sets are single-variant)
topvariant %>% slice_max(max_pp) %>% pull(variant_id) # Highest max PIP: chr1:77984833:C:A (max PIP = 1.00)
topvariant %>% slice_max(median_pp) %>% pull(variant_id) # Highest median PIP: chr1:77984833:C:A (median PIP = 0.962)

# For conditions without this variant, what is the lead?
non_coloc_cs = fm %>% 
    group_by(condition, credible_set) %>%
    summarise(coloc_var_in_cs = coloc_qtl_variant %in% variant_id) %>% 
    filter(!coloc_var_in_cs)

## A tibble: 2 × 3
## Groups:   condition [2]
#  condition credible_set coloc_var_in_cs
#  <chr>            <int> <lgl>          
#1 Plasma_ct            2 FALSE          
#2 T_blood              1 FALSE     


########
# Testing the DoE of the colocalising variant - subsetting to the conditions where this could be fine mapped
########
coloc_fm_sets = fm %>%
    filter(variant_id == coloc_qtl_variant) %>% 
    pull(condition)

nom_sumstats = read.delim(nomsumstat_file) %>% 
    filter(annotation %in% coloc_fm_sets)

sign_one = sign(nom_sumstats$slope[1])
all(sign(nom_sumstats$slope) == sign_one) # TRUE - so there are no opposite direction of effect. 