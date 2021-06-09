#!/bin/env Rscript

library(optparse)
library(Seurat)
library(stringr)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(data.table)
set.seed(10)

# path <- '/Volumes/My Passport/Gliomas/29.10.20/1_counts/Seq073_5/outs/filtered_feature_bc_matrix/'
# sample_name <- "Seq073_5"
# outpath <- '/Volumes/My Passport/Gliomas/15.02.21/2_getClusters/'
# min_genes <- 1000
# max_genes <- 7000
# perc_mito <- 10
# exclude <- NULL
# res <- 0.5
# normalization_method <- "CLR"

option_list = list(
  make_option(c("-i", "--inpath"), type="character", default=NULL,
              help="Cellranger output path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output path", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample name", metavar="character"),
  make_option(c("-r", "--resolution"), type="character", default=0.5,
              help="Resolution. Default 0.5", metavar="character"),
  make_option(c("-p", "--percMitochondria"), type="character", default=NULL,
              help="Maximum percentage of mitochondrial counts", metavar="character"),
  make_option(c("-m", "--minGenes"), type="numeric", default=500,
              help="Minimum number of genes detected per cell. Default 500.", metavar="numeric"),
  make_option(c("-M", "--maxGenes"), type="numeric", default=7000,
              help="Maximum number of genes detected per cell. Default 7000", metavar="numeric"),
  make_option(c("-n", "--normalizationMethod"), type="character", default="LogNormalize",
              help = "Seurat normalization method (LogNormalize | CLR)", metavar = "character"),
  make_option(c("-e", "--exclude"), type="character", default=NULL,
              help="Exclude cells", metavar="character"),
  make_option(c("-S", "--maxSize"), type="numeric", default=500,
              help = "Maximum size for global variables in MiB (future.globals.maxSize). This will increase your RAM usage.", metavar = "numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpath) | is.null(opt$outpath) | is.null(opt$sample)){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

options(future.globals.maxSize = as.numeric(as.character(opt$maxSize)) * 1024^2)
path <- ifelse(endsWith(opt$inpath, "/"), opt$inpath, paste(opt$inpath, '/', sep=''))
sample_name <- trimws(opt$sample)
outpath <- ifelse(endsWith(opt$outpath, "/"), opt$outpath, paste(opt$outpath, '/', sep=''))
min_genes <- opt$minGenes
max_genes <- opt$maxGenes
perc_mito <- opt$percMitochondria
exclude <- opt$exclude
res <- as.numeric(as.character(opt$resolution))
normalization_method <- as.character(opt$normalizationMethod)

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
print(paste("Resolution", as.character(res)))
print(paste("Normalization method", normalization_method))
print(paste("Future global max size", as.character(opt$maxSize)))

sample.data <- Read10X(data.dir = path)
sample <- CreateSeuratObject(counts = sample.data, project = sample_name, min.cells = 3, min.features = 200)
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
sample <- subset(sample, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes)
sample <- NormalizeData(sample, normalization.method = normalization_method, scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample)
sample <- RenameCells(sample, sample$orig.ident)

if(!is.null(exclude)){
  print(paste("sample used to have ", ncol(sample), "cells"))
  print(paste("We exclude ", length(excluding)))
  sample <- subset(sample, cells = colnames(sample)[which(!colnames(sample) %in% excluding)])
  print(paste("sample now has ", ncol(sample), "cells"))
}

if(!is.null(perc_mito)){
  sample <- subset(sample, subset = percent.mt < as.numeric(as.character(perc_mito)))
}

sample <- ScaleData(sample, features = all.genes)
sample <- RunPCA(sample, features = VariableFeatures(object = sample))
sample <- FindNeighbors(sample, dims = 1:10)
sample <- FindClusters(sample, resolution = res)
sample <- RunUMAP(sample, dims = 1:10)
write.csv(Embeddings(sample, reduction = "umap"), file = paste(outpath, '/', sample_name, "_cell_embeddings.csv", sep=''))

colours <- colorRampPalette(brewer.pal(8, "Accent"))(length(unique(sample$seurat_clusters)))
names(colours) <- as.character(unique(sample$seurat_clusters))
sample_colours <- merge(data.frame(seurat_clusters=unique(sample$seurat_clusters)), data.frame(cluster_colours = colours, seurat_clusters = names(colours)), by='seurat_clusters')
seurat_clusters_df <- as.data.frame(sample[["seurat_clusters"]])
seurat_clusters_df$Row.names <- rownames(seurat_clusters_df)
sample_colours <- merge(seurat_clusters_df, sample_colours, all.x = T)
rownames(sample_colours) <- sample_colours$Row.names
sample_colours <- structure(as.character(sample_colours$cluster_colours), names = as.character(sample_colours$Row.names))
sample <- AddMetaData(sample, sample_colours, col.name = "cluster_colours")

write.csv(sample[[c("seurat_clusters", "cluster_colours")]], file = paste(outpath, '/', sample_name, "_clusters.csv", sep=''))
sample[["cluster_colours"]] <- NULL

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

saveRDS(sample, file = paste(outpath, '/', sample_name, ".rds", sep=''))



