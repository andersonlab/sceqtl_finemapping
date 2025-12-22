suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(coloc))

##################
# Define options
##################
gene="ENSG00000143801"
sumstats.all.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/'
inputdir <- paste0("results/finemapping/", gene)
LDmatrix_path <- paste0("input/LDmatrix/", gene, ".ld")
LDsnps_path <- paste0("input/vcf/", gene, "_.z")

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

# Define celltypes
celltypes = do.call(rbind, cs) %>% pull(celltype) %>% unique()
names(cs) = celltypes

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
pairwise_coloc = function(nomlist, ct1, ct2){
    variants_1 = nomlist[[ct1]]$variant_id
    variants_2 = nomlist[[ct2]]$variant_id
    common = intersect(variants_1, variants_2)
    if(length(common) > 0){
        print(paste0(".. ", length(common), " variants in common between: ", ct1, " and ", ct2))
    } else {
        print(paste0("..No variants in common between: ", ct1, " and ", ct2))
        dummy=data.frame(
            nsnps = 0,
            hit1 = nomlist[[ct1]] %>% slice_min(pval_nominal) %>% pull(variant_id),
            hit2 = nomlist[[ct2]] %>% slice_min(pval_nominal) %>% pull(variant_id),
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
        LD = as.matrix(LDmatrix)   # Add LD matrix
    )

    D2 = list(
        beta = nomlist[[ct2]]$slope,
        varbeta = nomlist[[ct2]]$slope_se^2,
        snp = nomlist[[ct2]]$variant_id,
        position = as.numeric(unlist(strsplit(nomlist[[ct2]]$variant_id, "\\:"))[c(F,T,F,F)]),
        type = "quant",
        N = nomlist[[ct2]]$samplesize[1],  
        MAF = nomlist[[ct2]]$af,       
        LD = as.matrix(LDmatrix)   # Add LD matrix
    )

    res = coloc.susie(dataset1 = D1, dataset2 = D2)$summary
    res$cond1 = ct1
    res$cond2 = ct2
    rownames(res) = NULL
    return(res)
}

colocres = lapply(ctpairs, function(x){
    pairwise_coloc(nomlist=cs, ct1=x[1], ct2=x[2])
})
colocres = Filter(Negate(is.null), colocres) 
colocres = do.call(rbind, colocres) %>% 
    filter(!is.na(PP.H4.abf))



###################
# Run function to aggregate signals across pairwise comparisons and defines groups which are the same signal
###################