library(susieR)

args <- commandArgs(trailingOnly = TRUE)
gene=args[1]
nominalsumstats=args[2]
n=args[3]
LDmatrix_path=args[4]
LDsnps_path=args[5]
outpath=args[6]

#load the summary statistics for the relevant gene only
print("... Reading in summary statistics")
lines=readLines(nominalsumstats)
header=lines[1]
gene_lines=strsplit(grep(gene,lines,value=TRUE),'\t')
sumstats=as.data.frame(do.call(rbind,gene_lines))

# Make dummy - to be saved if finemapping is not possible
dummy <- data.frame(
        credible_set = character(),
        phenotype_id = character(),
        variant_id = character(),
        start_distance = numeric(),
        af = numeric(),
        ma_samples = integer(),
        ma_count = integer(),
        pval_nominal = numeric(),
        slope = numeric(),
        slope_se = numeric(),
        snp_number_in_window = integer(),
        PIP = numeric(),
        cs_log10bayesfactor = numeric(),
        cs_avg_r2 = numeric(),
        cs_min_r2 = numeric(),
        samplesize = integer(),
        stringsAsFactors = FALSE
    )   

if(nrow(sumstats) == 0){
    print("..Gene not found in sumstats, exiting and saving dummy output") 
    write.table(dummy, paste0(outpath,'crediblesets.txt'),row.names = FALSE, col.names = TRUE,quote=FALSE) #write results
    quit(save = "no", status = 0)
}
colnames(sumstats)=strsplit(header,'\t')[[1]]

#read in the ld matrix and .z file
LDmatrix=read.table(LDmatrix_path)
LDsnps=read.table(LDsnps_path,header=TRUE)

#Get index for LDmatrix SNPs not present in summary statistics
index=which(!(LDsnps$rsid %in% sumstats$variant_id))
print(paste(length(index), 'SNPs from LD matrix not present in summary statistics.'))

#only keep the SNPs from LDmatrix that are in the sumstats
print(paste('LD matrix cols/rows before filtering:',nrow(LDmatrix)))
print(paste('Rows of genotypes_window.z file before filtering:',nrow(LDsnps)))
if(length(index)>0){
    LDsnps_filt=LDsnps[-index,]
    LDmatrix_filt=LDmatrix[-index,-index]
} else {
    LDsnps_filt=LDsnps
    LDmatrix_filt=LDmatrix
}
rownames(LDsnps_filt) <- NULL #so the index is reset to 1,2,3,etc without gaps where the filtered out snps were
print(paste('LD matrix dimensions after filtering:',nrow(LDmatrix_filt)))
print(paste('Rows of genotypes_window.z filtered object:',nrow(LDsnps_filt)))



#Get index for summary statistic SNPs not present in LDmatrix
index_skipeQTLs=which(!(sumstats$variant_id %in% LDsnps$rsid))
skip_sumstats=sumstats[index_skipeQTLs,]
print(paste(length(index_skipeQTLs), 'SNPs from summary statistics not present in LD matrix.'))

#only keep sumstats snps present in LDmatrix
print(paste('SNPs in summary statistics before filtering:',nrow(sumstats)))
if(length(index_skipeQTLs)>0){
    sumstats_filt=sumstats[-index_skipeQTLs,]
    write.table(skip_sumstats,paste0(outpath,'missingSNPs.txt'),row.names = FALSE, col.names = TRUE, quote=FALSE)
    print("Outputting summary statistics for missing SNPs to missingSNPs.txt")
} else {
    sumstats_filt=sumstats
}
rownames(sumstats_filt) <- NULL #so the index is reset to 1,2,3,etc without gaps where the filtered out snps were
print(paste('SNPs in summary statistics after filtering:',nrow(sumstats_filt)))



#Remove any nas from sumstats and cases where se <= 0
nas=which(is.na(as.numeric(sumstats_filt$slope_se)) | is.na(as.numeric(sumstats_filt$slope)))
if(length(nas)>0){
    sumstats_filt=sumstats_filt[-nas,]
    LDmatrix_filt=LDmatrix_filt[-nas,-nas]
    LDsnps_filt=LDsnps_filt[-nas,]}

badses=which(as.numeric(sumstats_filt$slope_se) <= 0)
if(length(badses)>0){
    sumstats_filt=sumstats_filt[-badses,]
    LDmatrix_filt=LDmatrix_filt[-badses,-badses]
    LDsnps_filt=LDsnps_filt[-badses,]}



#check that the SNPs in sumstats and LDmatrix are the same and in the right order
if(identical(sumstats_filt$variant_id,LDsnps_filt$rsid)){
    print('SNPs in filtered summary stats and LD matrix match.')
} else {
    stop('SNPs in filtered summary stats and LD matrix DO NOT match.')
}

#make LDmatrix into an R matrix without colnames or rownames so it is symmetric for susie modelling
R=as.matrix(LDmatrix_filt)
rownames(R) <- NULL
colnames(R) <- NULL



#run susie
fitted_rss1 <- tryCatch({
    susie_rss(bhat = as.numeric(sumstats_filt$slope), #effect size
        shat = as.numeric(sumstats_filt$slope_se), #standard error
        n = as.numeric(n), #samplesize
        R = R, #LDmatrix – I assume of just the relevant SNPs?
        L = 10, #maximum number of causal variants, default=10
        estimate_residual_variance = FALSE) #TRUE because R is the in-sample LD matrix, but FALSE because this doesn't work for every gene
}, error = function(e) {
    print(paste("Error encountered in susie_rss:", e$message))
    print("Saving empty credible sets dataframe and exiting")
    write.table(dummy, paste0(outpath,'crediblesets.txt'),row.names = FALSE, col.names = TRUE,quote=FALSE)
    quit(save = "no", status = 0)
})



#prepare results
if(is.null(summary(fitted_rss1)$cs)){
    write.table(dummy, paste0(outpath,'crediblesets.txt'),row.names = FALSE, col.names = TRUE,quote=FALSE) #write results
} else {
    pips=summary(fitted_rss1)$vars #get posterior inclusion probabilities
    cs_pips=pips[pips$cs!=-1,] #keep only those from a credible set

    credible_sets=as.numeric(unlist(strsplit(summary(fitted_rss1)$cs$variable,','))) #get snp number of snps in credible sets
    cs_sumstats=sumstats_filt[credible_sets,] #get eQTL sumstats of snps in credible sets

        cs_sumstats_pips=merge(cs_sumstats,cs_pips,by='row.names') #merge sumstats and pip

    cs_susiesumstats=summary(fitted_rss1)$cs #get susie sumstats per credible set
    cs_susiesumstats$variable=NULL #remove useless column

    cs_finalsummary=merge(cs_sumstats_pips, cs_susiesumstats, by = "cs", all.x = TRUE) #get final summary of snps in credible sets with eQTL sumstats, PIPs and credible set sumstats
    cs_finalsummary$Row.names=NULL
    cs_finalsummary$samplesize=n
    colnames(cs_finalsummary)[colnames(cs_finalsummary) == 'variable']='snp_number_in_window' #rename columns
    colnames(cs_finalsummary)[colnames(cs_finalsummary) == 'variable_prob']='PIP'
    colnames(cs_finalsummary)[colnames(cs_finalsummary) == 'cs_log10bf']='cs_log10bayesfactor'
    colnames(cs_finalsummary)[colnames(cs_finalsummary) == 'cs']='credible_set'

#cs_log10bf = log10 bayes factor of comparing the solution of this model (cs independent credible sets) to cs -1 credible sets
#cs_avg_r2 = average r2 correlation between SNPs in the credible set
#cs_min_r2 = min r2 correlation between SNPs in the credible set

    write.table(cs_finalsummary,paste0(outpath,'crediblesets.txt'),row.names = FALSE, col.names = TRUE,quote=FALSE) #write results
}


