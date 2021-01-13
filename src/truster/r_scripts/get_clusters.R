#!/bin/env Rscript

library(optparse)
library(Seurat)
library(stringr)
library(dplyr)
library(patchwork)
set.seed(10)

option_list = list(
  make_option(c("-i", "--inpath"), type="character", default=NULL,
              help="Cellranger output path", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output path", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpath) | is.null(opt$outpath) | is.null(opt$sample)){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

path <- opt$inpath
sample_name <- opt$sample
outpath <- opt$outpath

#print(path)
#print(sample_name)

sample.data <- Read10X(data.dir = path)
sample <- CreateSeuratObject(counts = sample.data, project = sample_name, min.cells = 3, min.features = 200)
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
sample <- NormalizeData(sample, normalization.method = "LogNormalize", scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sample)
sample <- ScaleData(sample, features = all.genes)
sample <- RunPCA(sample, features = VariableFeatures(object = sample))
sample <- FindNeighbors(sample, dims = 1:10)
sample <- FindClusters(sample, resolution = 0.5)
sample <- RunUMAP(sample, dims = 1:10)

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

save(sample, file=paste(outpath, '/', sample_name, ".RData", sep=''))
