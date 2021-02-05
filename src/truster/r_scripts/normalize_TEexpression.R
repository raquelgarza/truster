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
# # 
mode = 'merged'
obj_name = "mergedCluster"
outdir = '/Volumes/My Passport/FetalCortex/01.02.21/3_mergeSamples/clusterPipeline/TEcountsNormalized/'
indir = '/Volumes/My Passport/FetalCortex/01.02.21/3_mergeSamples/clusterPipeline/TEcounts/'
rdata = '/Volumes/My Passport/FetalCortex/01.02.21/3_mergeSamples/fetalcortexPerSample.RData'

option_list = list(
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Merged samples or individual? (merged/individual)", metavar="character"),
  make_option(c("-r", "--RData"), type="character", default=NULL,
              help="Path to RData with Seurat object", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="Sample id or name", metavar="character"),
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="Path to TEcounts folder with output per cluster or parent directory from all samples folder", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Output directory for the matrix with normalized counts", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$indir) | is.null(opt$outdir) | is.null(opt$RData) | is.null(opt$mode) ){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

mode <- trimws(opt$mode)
rdata <- opt$RData
obj_name <- opt$name
indir <- ifelse(endsWith(opt$indir, "/"), opt$indir, paste(opt$indir, '/', sep=''))
outdir <- ifelse(endsWith(opt$outdir, "/"), opt$outdir, paste(opt$outdir, '/', sep=''))

if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

print(c("Mode: ", mode))
print(c("RData to be used: ", rdata))
print(c("Input path: ", indir))
print(c("Output path: ", outdir))

load(rdata)
if(mode == 'merged'){
  seurat.obj <- experiment
  cluster_sizes <- data.frame(cluster.size=table(seurat.obj@active.ident))  
  colnames(cluster_sizes) <- c('cluster', 'cluster.size')
  cluster_sizes$name <- paste(obj_name, cluster_sizes$cluster, sep='_')
}else{
  seurat.obj <- sample
  cluster_sizes <- data.frame(cluster.size=table(seurat.obj@active.ident))  
  colnames(cluster_sizes) <- c('cluster', 'cluster.size')
  cluster_sizes$name <- paste(obj_name, cluster_sizes$cluster, sep=".cluster_")
}

files <- list.files(indir, recursive = F)
coldata <- data.frame()
for(i in 1:length(files)){
  if(mode == 'merged'){
    cluster <- unlist(str_split(unlist(str_split(files[i], "mergedCluster_")[1])[2], "_"))[1]
    sample = paste(obj_name, cluster, sep="_")
    name <- paste(obj_name, cluster, sep="_")
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

num_reads <- data.frame(id=names(colSums(TEcounts[,rownames(coldata)])),
                        value=colSums(TEcounts[,rownames(coldata)]))
if(mode != "merged"){
  num_reads$sample_id <- sapply(str_split(num_reads$id, '[[.]]'), `[[`, 1)
  num_reads <- aggregate(num_reads$value, by=list(num_reads$sample_id), FUN=sum)
}

coldata <- merge(coldata, cluster_sizes[,c('name', 'cluster.size')], by='name')
rownames(coldata) <- coldata$name

colnames(num_reads) <- c('sample', 'num_reads')
rownames(num_reads) <- num_reads$sample

TEcounts <- subset(TEcounts, !startsWith(TEcounts$TE, 'ENSG'))

coldata <- merge(coldata, num_reads, by='sample')  
rownames(coldata) <- coldata$name
te_counts_size <- TEcounts[, rownames(coldata)]
te_counts_size_norm <- te_counts_size
te_counts_size_norm[] <- mapply('/', te_counts_size_norm[, rownames(coldata)], coldata$cluster.size)

if(mode != "merged"){
  te_counts_size_norm[] <- mapply('/', te_counts_size_norm[, rownames(coldata)], coldata$num_reads)  
  te_counts_size_norm <- te_counts_size_norm[, rownames(coldata)] * 1e+10
}

te_counts_size_norm$te_id <- rownames(te_counts_size_norm)
te_counts_size$te_id <- rownames(te_counts_size)
te_counts_size_norm_melt <- reshape2::melt(te_counts_size_norm, by=list(c('te_id')))
te_counts_size_melt <- reshape2::melt(te_counts_size, by=list(c('te_id')))

if(mode == "merged"){
  te_counts_size_norm_melt$cluster <- gsub("mergedCluster_", "", te_counts_size_norm_melt$variable)
  te_counts_size_melt$cluster <- gsub("mergedCluster_", "", te_counts_size_melt$variable)
}else{
  te_counts_size_norm_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_size_norm_melt$variable, '[[.]]'),`[[`, 2))
  te_counts_size_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_size_melt$variable, '[[.]]'),`[[`, 2))
}


te_counts_size_norm <- te_counts_size_norm[,c(colnames(te_counts_size_norm)[ncol(te_counts_size_norm)], colnames(te_counts_size_norm)[-ncol(te_counts_size_norm)])]
te_counts_size <- te_counts_size[,c(colnames(te_counts_size)[ncol(te_counts_size)], colnames(te_counts_size)[-ncol(te_counts_size)])]
fwrite(te_counts_size_norm, paste(outdir, 'TE_normalizedValues_matrix.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
# To compare with bulk data
fwrite(te_counts_size, paste(outdir, 'TE_rawValues_matrix.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)

if(mode == 'merged'){
  te_counts_size_norm_aggr <- aggregate(te_counts_size_norm_melt$value, list(te_counts_size_norm_melt$te_id, te_counts_size_norm_melt$variable, te_counts_size_norm_melt$cluster), mean)
  te_counts_size_aggr <- aggregate(te_counts_size_melt$value, list(te_counts_size_melt$te_id, te_counts_size_melt$variable, te_counts_size_melt$cluster), mean)
  # te_counts_size_aggr <- aggregate(te_counts_size_melt$value, list(te_counts_size_melt$te_id, te_counts_size_melt$cluster), mean)
  colnames(te_counts_size_norm_aggr) <- c('te_id', 'sample_cluster', 'cluster', 'value')
  colnames(te_counts_size_aggr) <- c('te_id', 'sample_cluster', 'cluster', 'value')
  fwrite(te_counts_size_norm_aggr, paste(outdir, 'TE_normalizedValues_aggregatedByClusters_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
  fwrite(te_counts_size_aggr, paste(outdir, 'TE_rawValues_aggregatedByClusters_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
}else{
  fwrite(te_counts_size_norm_melt, paste(outdir, 'TE_normalizedValues_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
  fwrite(te_counts_size_melt, paste(outdir, 'TE_rawValues_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
}



