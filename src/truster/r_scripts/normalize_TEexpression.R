#!/bin/env Rscript

library(optparse)
library(Seurat)
library(stringr)
library(data.table)
library(dplyr)
library(patchwork)
set.seed(10)

# Bed file of transposons
# Tab file of classification of transposons
# Path to RData 
# Path to TEcounts folder with output per cluster
# 
# 
# mode = 'per_sample'
# obj_name = 'Seq095_2'
# outdir = '/Volumes/My Passport/FetalCortex/06.05.21/2_getClusters/cluster_pipeline/TE_counts_normalized/multiple/Seq095_2/'
# indir = '/Volumes/My Passport/FetalCortex/06.05.21/2_getClusters/cluster_pipeline/TE_counts/multiple/Seq095_2/'
# rds = '/Volumes/My Passport/FetalCortex/06.05.21/2_getClusters/Seq095_2/Seq095_2.rds'
# # 
# mode = 'merged'
# obj_name = "gliomas"
# outdir = '/Volumes/My Passport/Gliomas/04.03.21/3_mergeSamples/clusterPipeline/TEcountsNormalized/multiple/'
# indir = '/Volumes/My Passport/Gliomas/04.03.21/3_mergeSamples/clusterPipeline/TEcounts/multiple/'
# rds = '/Volumes/My Passport/Gliomas/04.03.21/3_mergeSamples/gliomas.rds'

option_list = list(
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Merged samples or individual? (merged/individual)", metavar="character"),
  make_option(c("-r", "--RDS"), type="character", default=NULL,
              help="Path to RDS with Seurat object", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="Sample id or name", metavar="character"),
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="Path to TEcounts folder with output per cluster or parent directory from all samples folder", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Output directory for the matrix with normalized counts", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$indir) | is.null(opt$outdir) | is.null(opt$RDS) | is.null(opt$mode) ){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

mode <- trimws(opt$mode)
rds <- opt$RDS
obj_name <- opt$name
indir <- ifelse(endsWith(opt$indir, "/"), opt$indir, paste(opt$indir, '/', sep=''))
outdir <- ifelse(endsWith(opt$outdir, "/"), opt$outdir, paste(opt$outdir, '/', sep=''))

if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

print(c("Mode: ", mode))
print(c("RDS to be used: ", rds))
print(c("Input path: ", indir))
print(c("Output path: ", outdir))

seurat.obj <- readRDS(rds)

cluster_sizes <- data.frame(cluster.size=table(seurat.obj@active.ident))  
colnames(cluster_sizes) <- c('cluster', 'cluster.size')
cluster_sizes$name <- paste(obj_name, cluster_sizes$cluster, sep='.cluster_')

files <- list.files(indir, recursive = F)
coldata <- data.frame()
for(i in 1:length(files)){
  if(mode == 'merged'){
    # mergedCluster is how trusTEr is programmed to name the files
    cluster <- unlist(str_split(unlist(str_split(files[i], "mergedCluster_")[1])[2], "_"))[1]
    sample = paste(obj_name, cluster, sep="_")
    name <- paste(obj_name, cluster, sep=".cluster_")
    coldata <- rbind(coldata, data.frame(sample=sample, name = name, cluster=cluster))
  }else{
    file_name <- sub("_$", "", sapply(str_split(files[i], '.cntTable'), `[[`, 1))
    # sample <- unlist(str_split(file_name, '_clusters'))[1]
    cluster <- unlist(str_split(file_name, 'clusters_'))[2]
    name <- paste(obj_name, cluster, sep=".cluster_")
    coldata <- rbind(coldata, data.frame(sample=obj_name, name=name, cluster=cluster))
  }
  
  file <- paste(indir, files[i], sep='')
  
  if (i == 1){
    TEcounts <- fread(file, data.table = F)
    colnames(TEcounts) <- c("TE", name)
  }else{
    counts <- fread(file, data.table = F)
    colnames(counts) <- c("TE", name)
    TEcounts <- merge(TEcounts, counts, all=T, by="TE")
  }
}

rownames(TEcounts) <- TEcounts$TE
rownames(coldata) <- coldata$name

# num_reads <- data.frame(id=names(colSums(TEcounts[which(startsWith(TEcounts$TE, "ENS")),rownames(coldata)])),
#                         value=colSums(TEcounts[which(startsWith(TEcounts$TE, "ENS")),rownames(coldata)]))
# if(mode != "merged"){
#   num_reads$sample_id <- sapply(str_split(num_reads$id, '[[.]]'), `[[`, 1)
#   # If we are normalizing by seq depth per sample, we can aggregate per sample id
#   num_reads <- aggregate(num_reads$value, by=list(num_reads$sample_id), FUN=sum)
# }

coldata <- merge(coldata, cluster_sizes[,c('name', 'cluster.size')], by='name')
rownames(coldata) <- coldata$name

cluster_total_num_reads <- function(cluster){
  return(sum(rowSums(FetchData(subset(seurat.obj, seurat_clusters == cluster), vars=rownames(seurat.obj)))))
}
coldata$num_reads <- sapply(coldata$cluster, FUN = cluster_total_num_reads)

# colnames(num_reads) <- c('sample', 'num_reads')
# rownames(num_reads) <- num_reads$sample

TEcounts <- subset(TEcounts, !startsWith(TEcounts$TE, 'ENS'))

# if(mode != "merged"){
#   coldata <- merge(coldata, num_reads, by='sample')  
# }

rownames(coldata) <- coldata$name
te_counts <- TEcounts[, rownames(coldata)]
te_counts_size_norm <- TEcounts[, rownames(coldata)]
te_counts_size_norm[] <- mapply('/', te_counts_size_norm[, rownames(coldata)], coldata$cluster.size)

# if(mode != "merged"){
#   # ONLY if it's per sample, divide the reads by the sequencing depth of the sample 
#   # Otherwise we end up dividing per number of reads in a cluster which is against what we want...
#   # Please ensure to take away "num_reads" if we are in merged mode
#   te_counts_size_norm[] <- mapply('/', te_counts_size_norm[, rownames(coldata)], coldata$num_reads)  
#   te_counts_size_norm <- te_counts_size_norm[, rownames(coldata)] * 1e+10
# }

te_counts_size_norm$te_id <- rownames(te_counts_size_norm)
te_counts_size_norm_melt <- reshape2::melt(te_counts_size_norm, by=list(c('te_id')))
te_counts_size_norm_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_size_norm_melt$variable, '[[.]]'),`[[`, 2))

cells_clusters <- data.frame(seurat.obj$seurat_clusters)
colnames(cells_clusters) <- "cluster"
cells_clusters$cell_ids <- rownames(cells_clusters)
te_counts_size_norm_melt_percell <- merge(cells_clusters, te_counts_size_norm_melt, by='cluster', all = T)
te_counts_size_norm_melt_percell <- reshape2::dcast(te_counts_size_norm_melt_percell, formula = te_id~cell_ids, value=value)
rownames(te_counts_size_norm_melt_percell) <- te_counts_size_norm_melt_percell$te_id  
te_counts_size_norm_melt_percell <- te_counts_size_norm_melt_percell[,-1]

# Norm by cluster size
seurat.obj[["TE_norm_cluster_size"]] <- CreateAssayObject(counts = te_counts_size_norm_melt_percell)

te_counts$te_id <- rownames(te_counts)
te_counts_melt <- reshape2::melt(te_counts, by=list(c('te_id')))
te_counts_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_melt$variable, '[[.]]'),`[[`, 2))
te_counts_melt_percell <- merge(cells_clusters, te_counts_melt, by='cluster', all = T)
te_counts_melt_percell <- reshape2::dcast(te_counts_melt_percell, formula = te_id~cell_ids, value=value)
rownames(te_counts_melt_percell) <- te_counts_melt_percell$te_id  
te_counts_melt_percell <- te_counts_melt_percell[,-1]

seurat.obj[["TE_raw"]] <- CreateAssayObject(counts = te_counts_melt_percell)

te_counts_size_reads_norm <- te_counts_size_norm[, rownames(coldata)]
te_counts_size_reads_norm[] <- mapply('/', te_counts_size_norm[, rownames(coldata)], coldata$num_reads)
te_counts_size_reads_norm <- te_counts_size_reads_norm[, rownames(coldata)] * 1e+7

te_counts_size_reads_norm$te_id <- rownames(te_counts_size_reads_norm)
te_counts_size_reads_norm_melt <- reshape2::melt(te_counts_size_reads_norm, by=list(c('te_id')))
te_counts_size_reads_norm_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_size_reads_norm_melt$variable, '[[.]]'),`[[`, 2))

te_counts_size_reads_norm_melt_percell <- merge(cells_clusters, te_counts_size_reads_norm_melt, by='cluster', all = T)
te_counts_size_reads_norm_melt_percell <- reshape2::dcast(te_counts_size_reads_norm_melt_percell, formula = te_id~cell_ids, value=value)
rownames(te_counts_size_reads_norm_melt_percell) <- te_counts_size_reads_norm_melt_percell$te_id  
te_counts_size_reads_norm_melt_percell <- te_counts_size_reads_norm_melt_percell[,-1]

# Norm by cluster size and sequencing depth
seurat.obj[["TE_norm_cluster_size_num_reads"]] <- CreateAssayObject(counts = te_counts_size_reads_norm_melt_percell)

te_counts_size_norm <- te_counts_size_norm[,c(colnames(te_counts_size_norm)[ncol(te_counts_size_norm)], colnames(te_counts_size_norm)[-ncol(te_counts_size_norm)])]
te_counts_size_reads_norm <- te_counts_size_reads_norm[,c(colnames(te_counts_size_reads_norm)[ncol(te_counts_size_reads_norm)], colnames(te_counts_size_reads_norm)[-ncol(te_counts_size_reads_norm)])]
te_counts <- te_counts[,c(colnames(te_counts)[ncol(te_counts)], colnames(te_counts)[-ncol(te_counts)])]
fwrite(te_counts_size_norm, paste(outdir, 'TE_norm_cluster_size_matrix.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
fwrite(te_counts_size_reads_norm, paste(outdir, 'TE_norm_cluster_size_num_reads_matrix.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
# To compare with bulk data
fwrite(te_counts, paste(outdir, 'TE_raw_matrix.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)

# if(mode == 'merged'){
#   te_counts_size_norm_aggr <- aggregate(te_counts_size_norm_melt$value, list(te_counts_size_norm_melt$te_id, te_counts_size_norm_melt$variable, te_counts_size_norm_melt$cluster), mean)
#   te_counts_size_aggr <- aggregate(te_counts_size_melt$value, list(te_counts_size_melt$te_id, te_counts_size_melt$variable, te_counts_size_melt$cluster), mean)
#   colnames(te_counts_size_norm_aggr) <- c('te_id', 'sample_cluster', 'cluster', 'value')
#   colnames(te_counts_size_aggr) <- c('te_id', 'sample_cluster', 'cluster', 'value')
#   fwrite(te_counts_size_norm_aggr, paste(outdir, 'TE_normalizedValues_aggregatedByClusters_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
#   fwrite(te_counts_size_aggr, paste(outdir, 'TE_rawValues_aggregatedByClusters_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
# }else{
#   fwrite(te_counts_size_norm_melt, paste(outdir, 'TE_normalizedValues_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
#   fwrite(te_counts_size_melt, paste(outdir, 'TE_rawValues_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
# }


saveRDS(seurat.obj, file = rds)



