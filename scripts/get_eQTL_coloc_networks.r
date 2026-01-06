suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(digest))

##################
# Define functions
##################
pairwise_coloc = function(nomlist, ct1, ct2){
    variants_1 = nomlist[[ct1]]$variant_id
    variants_2 = nomlist[[ct2]]$variant_id
    common = intersect(variants_1, variants_2)
    if(length(common) == 0){
        dummy=data.frame(
            nsnps = 0,
            hit1 = nomlist[[ct1]] %>% slice_min(pval_nominal, with_ties = FALSE) %>% pull(variant_id),
            hit2 = nomlist[[ct2]] %>% slice_min(pval_nominal, with_ties = FALSE) %>% pull(variant_id),
            PP.H0.abf = 0,
            PP.H1.abf = 0,
            PP.H2.abf = 0,
            PP.H3.abf = 1,
            PP.H4.abf = 0,
            idx1 = 1,
            idx2 = 1, 
            cond1 = ct1, 
            cond2 = ct2
        )
        return(dummy)
    }

    D1 = list(
        beta = nomlist[[ct1]]$slope,
        varbeta = nomlist[[ct1]]$slope_se^2,
        snp = nomlist[[ct1]]$variant_id,
        position = as.numeric(unlist(strsplit(nomlist[[ct1]]$variant_id, "\\:"))[c(F,T,F,F)]),
        type = "quant",
        N = nomlist[[ct1]]$samplesize[1],  
        MAF = nomlist[[ct1]]$af,       
        LD = as.matrix(LDmatrix)   
    )

    D2 = list(
        beta = nomlist[[ct2]]$slope,
        varbeta = nomlist[[ct2]]$slope_se^2,
        snp = nomlist[[ct2]]$variant_id,
        position = as.numeric(unlist(strsplit(nomlist[[ct2]]$variant_id, "\\:"))[c(F,T,F,F)]),
        type = "quant",
        N = nomlist[[ct2]]$samplesize[1],  
        MAF = nomlist[[ct2]]$af,       
        LD = as.matrix(LDmatrix)   
    )

    res = suppressMessages(coloc.susie(dataset1 = D1, dataset2 = D2)$summary)
    res$cond1 = ct1
    res$cond2 = ct2
    rownames(res) = NULL
    return(res)
}

get_cs_summary <- function(cs_list) {
  cs_summary <- do.call(rbind, cs_list) %>%
    group_by(credible_set, phenotype_id, celltype) %>%
    summarise(
      size = n(),
      maxpip = max(PIP, na.rm = TRUE),
      maxpip_variant = variant_id[which.max(PIP)],
      .groups = "drop"
    ) %>% 
    rename(
        idx = credible_set
    ) %>%
    mutate(
        hash = substr(sapply(paste(idx, phenotype_id, celltype, sep = "_"), function(x) digest(x, algo = "md5")), 1, 10)
    )
    
  return(cs_summary)
}

merge_shared_effects <- function(cs_summary, coloc_res, colocthresh, method = "soft"){
    coloc_filtered = colocres %>%
        filter(PP.H4.abf > colocthresh)

    if(nrow(coloc_filtered) > 0){
        print("..Colocalised credible sets found, merging")
        # Create edge list for network
        edges = coloc_filtered %>%
            select(hash1, hash2, PP.H4.abf, nsnps) %>%
            rename(from = hash1, to = hash2, weight = PP.H4.abf)
        
        # Create node attributes
        nodes = bind_rows(
            coloc_filtered %>% select(hash1, cond1, idx1) %>%
                rename(hash = hash1, celltype = cond1, credible_set = idx1),
            coloc_filtered %>% select(hash2, cond2, idx2) %>%
                rename(hash = hash2, celltype = cond2, credible_set = idx2)
        ) %>%
            distinct()
        
        # Create igraph object
        g = graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
            
        # Identify connected components (unique colocalization sets)
        comps = components(g)
        node_components = nodes %>%
            rename(name = hash) %>%
            left_join(data.frame(name = V(g)$name, component = comps$membership), by = "name")
        coloc_sets = split(node_components, node_components$component)
        print(paste0("..Found ", length(coloc_sets), " unique colocalization sets (connected components)"))

        if(method == "hard"){
            # For hard method, only merge if all nodes in component are fully connected (complete subgraph/clique)
            valid_components = sapply(coloc_sets, function(comp_nodes){
                comp_hashes = comp_nodes$name
                if(length(comp_hashes) == 1) return(TRUE)
                
                # Get subgraph for this component
                subg = induced_subgraph(g, comp_hashes)
                
                # Check if it's a complete graph (all nodes connected to all other nodes)
                n_nodes = vcount(subg)
                n_edges = ecount(subg)
                expected_edges = n_nodes * (n_nodes - 1) / 2
                
                return(n_edges == expected_edges)
            })
            
            # Filter to only valid components
            valid_comp_ids = as.numeric(names(coloc_sets)[valid_components])
            node_components = node_components %>%
                mutate(component = ifelse(component %in% valid_comp_ids, component, NA))
            
            print(paste0("..Found ", sum(valid_components), " fully connected components (hard method)"))
        }

        # Create deterministic component labels from the sorted hashes within each component
        valid_node_components = node_components[!is.na(node_components$component), ]
        
        if(nrow(valid_node_components) > 0){
            comp_labels = sapply(split(valid_node_components$name, valid_node_components$component), function(x){
                paste(sort(unique(x)), collapse = "|")
            })
            
            # Long mapping: each hash -> its collapsed component label
            comp_labels_df = data.frame(
                hash = unlist(strsplit(comp_labels, "\\|")),
                collapsed_hash = rep(unname(comp_labels), times = sapply(strsplit(comp_labels, "\\|"), length)),
                stringsAsFactors = FALSE
            )
        } else {
            # No valid components, create empty dataframe
            comp_labels_df = data.frame(
                hash = character(),
                collapsed_hash = character(),
                stringsAsFactors = FALSE
            )
        }
        
        # Merge with summary and flatten
        hash_col = ifelse(method == "hard", "strict_collapsed_hash", "collapsed_hash")
        cs_summary = cs_summary %>%
            left_join(comp_labels_df, by = c("hash" = "hash")) %>%
            rename(!!hash_col := collapsed_hash) %>%
            mutate(
                !!hash_col := ifelse(is.na(.data[[hash_col]]), hash, .data[[hash_col]])
            )
        
    } else {
        print("..No colocalizations with PP.H4.abf > 0.75 found")
        hash_col = ifelse(method == "hard", "strict_collapsed_hash", "collapsed_hash")
        cs_summary[[hash_col]] = cs_summary$hash
    }

    return(cs_summary)
}

##################
# Define options
##################
args <- commandArgs(trailingOnly = TRUE)
gene = args[1] # gene="ENSG00000143801" # PSEN2
inputdir <- args[2] #inputdir <- paste0("results/finemapping/", gene)
LDmatrix_path = args[3] # LDmatrix_path <- paste0("input/LDmatrix/", gene, ".ld")
LDsnps_path = args[4] # LDsnps_path <- paste0("input/vcf/", gene, "_.z")
colocthresh = as.numeric(args[5]) # colocthresh <- 0.75

###################
# Define the celltypes we are going to run coloc between, load the credible sets
###################
print("..Loading credible sets")
inputf = list.files(inputdir, pattern=".crediblesets", full.names=T)
cs = lapply(inputf, function(x){
    temp = read.delim(x, sep = " ")
    if(nrow(temp) > 1){
        temp$celltype = strsplit(x, "\\/") %>% unlist() %>% tail(1) %>% gsub(".crediblesets.txt", "", .)
        return(temp)
    }
})
cs = Filter(Negate(is.null), cs) 
print(paste0("..Found credible sets in ", length(cs), " conditions"))

if(length(cs) == 0){
    print("..No credible sets found, saving dummy output and quitting")
    dummy = data.frame(
        credible_set = character(),
        phenotype_id = character(),
        variant_id = character(),
        start_distance = numeric(),
        af = numeric(),
        ma_samples = numeric(),
        ma_count = numeric(),
        pval_nominal = numeric(),
        slope = numeric(),
        slope_se = numeric(),
        snp_number_in_window = integer(),
        PIP = numeric(),
        cs_log10bayesfactor = numeric(),
        cs_avg_r2 = numeric(),
        cs_min_r2 = numeric(),
        samplesize = integer(),
        celltype = character(),
        size = integer(),
        maxpip = numeric(),
        maxpip_variant = character(),
        hash = character(),
        collapsed_hash = character(),
        strict_collapsed_hash = character(),
        stringsAsFactors = FALSE
        )
    write.table(dummy, file = paste0("results/coloc/", gene, "_coloc_postfinemap.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

    # Empty pairwise coloc results
    dummy_coloc = data.frame(
        nsnps = integer(),
        hit1 = character(),
        hit2 = character(),
        PP.H0.abf = numeric(),
        PP.H1.abf = numeric(),
        PP.H2.abf = numeric(),
        PP.H3.abf = numeric(),
        PP.H4.abf = numeric(),
        idx1 = integer(),
        idx2 = integer(),
        cond1 = character(),
        cond2 = character(),
        hash1 = character(),
        hash2 = character(),
        stringsAsFactors = FALSE
    )
    print("..Also saving dummy coloc")
    write.table(dummy_coloc, file = paste0("results/coloc/", gene, "_pairwise_coloc_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    quit(save = "no", status = 0)
}

# Define celltypes
celltypes = do.call(rbind, cs) %>% pull(celltype) %>% unique()
names(cs) = celltypes

if(length(cs) < 2){
    print("..Saving dummy, nothing to merge")
    if(length(cs) == 1){
        dummy = do.call(rbind, cs)
        dummy$size = nrow(dummy)
        dummy$maxpip = max(dummy$PIP)
        dummy$maxpip_variant = dummy %>% slice_max(PIP, with_ties = FALSE) %>% pull(variant_id)
        dummy$hash = "a1"
        dummy$collapsed_hash = "a1"
        dummy$strict_collapsed_hash = "a1"
        
    } else {
        dummy = data.frame(
            credible_set = character(),
            phenotype_id = character(),
            variant_id = character(),
            start_distance = numeric(),
            af = numeric(),
            ma_samples = numeric(),
            ma_count = numeric(),
            pval_nominal = numeric(),
            slope = numeric(),
            slope_se = numeric(),
            snp_number_in_window = integer(),
            PIP = numeric(),
            cs_log10bayesfactor = numeric(),
            cs_avg_r2 = numeric(),
            cs_min_r2 = numeric(),
            samplesize = integer(),
            celltype = character(),
            size = integer(),
            maxpip = numeric(),
            maxpip_variant = character(),
            hash = character(),
            collapsed_hash = character(),
            strict_collapsed_hash = character(),
            stringsAsFactors = FALSE
            )
    }
    write.table(dummy, file = paste0("results/coloc/", gene, "_coloc_postfinemap.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

    # Empty pairwise coloc results
    dummy_coloc = data.frame(
        nsnps = integer(),
        hit1 = character(),
        hit2 = character(),
        PP.H0.abf = numeric(),
        PP.H1.abf = numeric(),
        PP.H2.abf = numeric(),
        PP.H3.abf = numeric(),
        PP.H4.abf = numeric(),
        idx1 = integer(),
        idx2 = integer(),
        cond1 = character(),
        cond2 = character(),
        hash1 = character(),
        hash2 = character(),
        stringsAsFactors = FALSE
    )
    print("..Also saving dummy coloc")
    write.table(dummy_coloc, file = paste0("results/coloc/", gene, "_pairwise_coloc_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    quit(save = "no", status = 0)
}
# Define comparisons
ctpairs = combn(celltypes, 2, simplify = FALSE)

###################
# Get LD
###################
print("..Loading LD")
LDmatrix = read.table(LDmatrix_path)
LDsnps = read.table(LDsnps_path,header=TRUE)
dimnames(LDmatrix) = list(LDsnps$rsid, LDsnps$rsid)

###################
# Run coloc 
###################
print("..Running pairwise coloc between celltypes")
total_pairs = length(ctpairs)
progress_interval = max(1, ceiling(total_pairs * 0.05))
start_time = Sys.time()
last_checkpoint_time = start_time

colocres = lapply(seq_along(ctpairs), function(i){
    # Run coloc
    result = pairwise_coloc(nomlist=cs, ct1=ctpairs[[i]][1], ct2=ctpairs[[i]][2])
    
    # Progress tracking every 5%
    if (i %% progress_interval == 0 || i == total_pairs) {
        current_time = Sys.time()
        elapsed_total = difftime(current_time, start_time, units = "secs")
        elapsed_last = difftime(current_time, last_checkpoint_time, units = "secs")
        percent_done = round((i / total_pairs) * 100, 1)
        print(paste0("....Progress: ", percent_done, "% (", i, "/", total_pairs, 
                     "), Time for last checkpoint: ", round(elapsed_last, 2), "s, Total time: ", 
                     round(elapsed_total, 2), "s"))
        last_checkpoint_time <<- current_time
    }
    
    return(result)
})

colocres = Filter(Negate(is.null), colocres) 
colocres = bind_rows(colocres) %>% 
    filter(!is.na(PP.H4.abf))

# Add hash
colocres = colocres %>%
    mutate(
        hash1 = substr(sapply(paste(idx1, sym(gene), cond1, sep = "_"), function(x) digest(x, algo = "md5")), 1, 10),
        hash2 = substr(sapply(paste(idx2, sym(gene), cond2, sep = "_"), function(x) digest(x, algo = "md5")), 1, 10)
    )   

###################
# Use this to aggregate effects and define independent eQTL signals
###################
cs_summary = get_cs_summary(cs) # hashes will match those of above

###################
# Create network from colocalization results, merge with credible set info
###################
print("..Merging colocalised credible sets - SOFT")
shared_cs = merge_shared_effects(cs_summary, colocres, colocthresh, method = "soft")
print("..Merging colocalised credible sets - HARD")
shared_cs_hard = merge_shared_effects(cs_summary, colocres, colocthresh, method = "hard")

cs_all = bind_rows(cs) %>%
    left_join(shared_cs, by = c("credible_set" = "idx", "phenotype_id", "celltype")) %>% 
    left_join(shared_cs_hard %>% select(idx, phenotype_id, celltype, strict_collapsed_hash), 
              by = c("credible_set" = "idx", "phenotype_id", "celltype"))

ncelltypes = length(unique(cs_all$celltype))
print(paste0("-- From ", ncelltypes, " celltypes:"))

nhash = length(unique(cs_all$collapsed_hash))
print(paste0("..Have found ", nhash, " unique sets after merging (soft)"))

nhash_hard = length(unique(cs_all$strict_collapsed_hash))
print(paste0("..Have found ", nhash_hard, " unique sets after merging (hard)"))

###################
# Save
###################
print("..Saving results")
write.table(cs_all, file = paste0("results/coloc/", gene, "_coloc_postfinemap.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(colocres, file = paste0("results/coloc/", gene, "_pairwise_coloc_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
print("..DONE!")



