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

# Rscript r_scripts/normalize_TEexpression.R -m merged -r /Volumes/LaCie/Gliomas/01.01.21/3_mergeSamplesgliomas.RData -i /Volumes/LaCie/Gliomas/01.01.21/3_mergeSamples/clusterPipeline/TEcounts/ -o /Volumes/LaCie/Gliomas/01.01.21/3_mergeSamples/clusterPipeline/TEcountsNorm/
# rdata <- '/Volumes/LaCie/Gliomas/01.01.21/3_mergeSamplesgliomas.RData'
# indir <- '/Volumes/LaCie/Gliomas/01.01.21/3_mergeSamples/clusterPipeline/TEcounts/'
# rdata <- "/Volumes/LaCie/FetalCortex/2_getClusters/DA094/DA094.RData"
# indir <- "/Volumes/LaCie/FetalCortex/2_getClusters/clusterPipeline/TEcounts/DA094/"

option_list = list(
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Merged samples or individual? (merged/individual)", metavar="character"),
  make_option(c("-r", "--RData"), type="character", default=NULL,
              help="Path to RData with Seurat object", metavar="character"),
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

mode <- opt$mode
rdata <- opt$RData
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
  cluster_sizes <- data.frame(cluster_size=table(experiment@active.ident))  
}else{
  cluster_sizes <- data.frame(cluster_size=table(sample@active.ident))
}
colnames(cluster_sizes) <- c('cluster', 'cluster_size')

files <- list.files(indir, recursive = T)
coldata <- data.frame()
for(i in 1:length(files)){
  if(mode == 'merged'){
    sample <- unlist(str_split(files[i], '/'))[1]
    cluster <- unlist(str_split(unlist(str_split(unlist(str_split(files[i], '/'))[2], "clusters_")[1])[2], "_"))[1]
  }else{
    sample <- unlist(str_split(files[i], '_'))[2]
    cluster <- unlist(str_split(files[i], '_'))[4]
  }
  coldata <- rbind(coldata, data.frame(sample=sample, cluster=cluster))
  name <- paste(sample, cluster, sep=".cluster_")
  
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

coldata$name <- paste(coldata$sample, coldata$cluster, sep=".cluster_")
cluster_sizes <- data.frame(cluster_size=table(experiment@active.ident))
colnames(cluster_sizes) <- c('cluster', 'cluster_size')

coldata <- merge(coldata, cluster_sizes, by='cluster')
rownames(coldata) <- coldata$name

num_reads <- data.frame(value=colSums(TEcounts[,rownames(coldata)]), id=names(colSums(TEcounts[,rownames(coldata)])))
num_reads$sample_id <- sapply(str_split(num_reads$id, '[[.]]'), `[[`, 1)
num_reads <- aggregate(num_reads$value, by=list(num_reads$sample_id), FUN=sum)
colnames(num_reads) <- c('sample', 'num_reads')
rownames(num_reads) <- num_reads$sample

TEcounts <- subset(TEcounts, !startsWith(TEcounts$TE, 'ENSG'))

coldata <- merge(coldata, num_reads, by='sample')
rownames(coldata) <- coldata$name
te_counts_size <- TEcounts[, rownames(coldata)]
te_counts_size[] <- mapply('/', te_counts_size[, rownames(coldata)], coldata$cluster_size)
te_counts_size[] <- mapply('/', te_counts_size[, rownames(coldata)], coldata$num_reads)
te_counts_size <- te_counts_size[, rownames(coldata)] * 1e+10

te_counts_size$te_id <- rownames(te_counts_size)
te_counts_size_melt <- reshape2::melt(te_counts_size, by=list(c('te_id')))
te_counts_size_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_size_melt$variable, '[[.]]'),`[[`, 2))

te_counts_size <- te_counts_size[,c(colnames(te_counts_size)[ncol(te_counts_size)], colnames(te_counts_size)[-ncol(te_counts_size)])]
fwrite(te_counts_size, paste(outdir, 'TE_normalizedValues_matrix.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)

if(mode == 'merged'){
  te_counts_size_aggr <- aggregate(te_counts_size_melt$value, list(te_counts_size_melt$te_id, te_counts_size_melt$cluster), mean)
  colnames(te_counts_size_aggr) <- c('te_id', 'cluster', 'value')
  fwrite(te_counts_size_aggr, paste(outdir, 'TE_normalizedValues_aggregatedByClusters_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
}else{
  fwrite(te_counts_size_melt, paste(outdir, 'TE_normalizedValues_melted.csv', sep=''), quote = F, row.names = F, col.names = T, verbose = T)
}
