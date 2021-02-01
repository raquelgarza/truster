#!/bin/env Rscript

library(optparse)
library(Seurat)
library(stringr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(data.table)
set.seed(10)

option_list = list(
  make_option(c("-i", "--inpath"), type="character", default=NULL,
              help="Cellranger output path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output path", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample name", metavar="character"),
  make_option(c("-p", "--percMitochondria"), type="character", default=NULL,
              help="Maximum percentage of mitochondrial counts", metavar="character"),
  make_option(c("-m", "--minGenes"), type="character", default=500,
              help="Minimum number of genes detected per cell. Default 500.", metavar="character"),
  make_option(c("-e", "--exclude"), type="character", default=NULL,
              help="Exclude cells", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpath) | is.null(opt$outpath) | is.null(opt$sample)){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

path <- ifelse(endsWith(opt$inpath, "/"), opt$inpath, paste(opt$inpath, '/', sep=''))
sample_name <- trimws(opt$sample)
outpath <- ifelse(endsWith(opt$outpath, "/"), opt$outpath, paste(opt$outpath, '/', sep=''))
min_genes <- opt$minGenes
perc_mito <- opt$percMitochondria
exclude <- opt$exclude

if(!is.null(exclude)){
  excluding <- fread(exclude, data.table = F, header = F)[,1]
  print(paste("Excluding some cells from ", exclude))
}
print(paste("Path", path))
print(paste("Sample", sample_name))
print(paste("Outpath", outpath))
print(paste("Min genes", min_genes))
print(paste("Perc mito", perc_mito))
print(paste("Exclude", exclude))

sample.data <- Read10X(data.dir = path)
sample <- CreateSeuratObject(counts = sample.data, project = sample_name, min.cells = 3, min.features = 200)
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample)

if(!is.null(exclude)){
  print(paste("sample used to have ", ncol(sample), "cells"))
  print(paste("We exclude ", length(excluding)))
  sample <- subset(sample, cells = colnames(sample)[which(!colnames(sample) %in% excluding)])
  print(paste("sample now has ", ncol(sample), "cells"))
}

if(!is.null(perc_mito)){
  sample <- subset(sample, subset = percent.mt < as.numeric(as.character(perc_mito)))
}

sample <- subset(sample, subset = nCount_RNA > as.numeric(as.character(min_genes)))
sample <- ScaleData(sample, features = all.genes)
sample <- RunPCA(sample, features = VariableFeatures(object = sample))
sample <- FindNeighbors(sample, dims = 1:10)
sample <- FindClusters(sample, resolution = 0.5)
sample <- RunUMAP(sample, dims = 1:10)
write.csv(Embeddings(sample, reduction = "umap"), file = paste(outpath, '/', sample_name, "_cell_embeddings.csv", sep=''))

colours <- colorRampPalette(brewer.pal(8, "Accent"))(length(unique(sample$seurat_clusters)))
names(colours) <- as.character(unique(sample$seurat_clusters))
sample@meta.data$cellIds <- rownames(sample@meta.data)
sample_colours <- merge(data.frame(seurat_clusters=unique(sample$seurat_clusters)), data.frame(cluster_colours = colours, seurat_clusters = names(colours)), by='seurat_clusters')
sample@meta.data <- merge(sample@meta.data, sample_colours, by='seurat_clusters')
rownames(sample@meta.data) <- sample@meta.data$cellIds
write.csv(sample@meta.data[,c('seurat_clusters', 'cluster_colours'), drop=F], file = paste(outpath, '/', sample_name, "_clusters.csv", sep=''))

df <- as.data.frame(sample$seurat_clusters)
colnames(df) <- 'clusters'

for (k in 1:length(unique(df$clusters))){
    cluster <- unique(df$clusters)[k]
    df.cluster <- subset(df, df$clusters == cluster)
    df.cluster$barcode <- rownames(df.cluster)

    file_name <- paste(outpath, '/', sample_name, '_clusters_', cluster, '.tsv', sep = '')
    print(df.cluster$barcode)
    write.table(df.cluster$barcode, file = file_name, row.names = F, col.names = F, quote = F)
}

if(!is.null(exclude)){
  save(sample, excluding, file = paste(outpath, '/', sample_name, ".RData", sep=''))
}else{
  save(sample, file = paste(outpath, '/', sample_name, ".RData", sep=''))  
}



