#!/bin/env Rscript

library(optparse)
library(Seurat)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(patchwork)
set.seed(10)

# paths <- c('/Volumes/My Passport/FetalCortex/16.01.21/2_getClustersExclusive/DA094/DA094.RData',
#            '/Volumes/My Passport/FetalCortex/16.01.21/2_getClustersExclusive/DA103/DA103.RData',
#            '/Volumes/My Passport/FetalCortex/16.01.21/2_getClustersExclusive/DA140/DA140.RData',
#            '/Volumes/My Passport/FetalCortex/16.01.21/2_getClustersExclusive/Seq098_2/Seq098_2.RData')
# 
# ids <- c("DA094", "DA103", "DA140", "Seq098_2")
# outpath <- "/Volumes/My Passport/FetalCortex/16.01.21/3_mergeSamples/"
# experiment_name <- "fetalcortex"

option_list = list(
  make_option(c("-i", "--inpath"), type="character", default=NULL,
              help="RData paths", metavar="character"),
  make_option(c("-s", "--ids"), type="character", default=NULL,
              help="Sample ids in order of appearance in -i", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL,
              help="Output path", metavar="character"),
  make_option(c("-e", "--experimentName"), type="character", default=NULL,
              help="Experiment name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpath) | is.null(opt$outpath) ){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

paths <- trimws(unlist(str_split(opt$inpath, ',')))
ids <- trimws(unlist(str_split(opt$ids, ',')))
outpath <- ifelse(endsWith(opt$outpath, "/"), opt$outpath, paste(opt$outpath, '/', sep=''))
experiment_name <- trimws(opt$experimentName)

print(c("Input paths: ", paths))
print(c("Input ids: ", ids))
print(c("Output path: ", outpath))
print(c("Experiment name: ", experiment_name))

for(i in 1:length(paths)){
  path <- paths[i]

  load(path)

  if(i == 1)
  {
    experiment <- sample
  }
  else
  {
    experiment <- merge(experiment, sample)
  }
}

experiment[["percent.mt"]] <- PercentageFeatureSet(experiment, pattern = "^MT-")
experiment <- NormalizeData(experiment, normalization.method = "LogNormalize", scale.factor = 10000)
experiment <- FindVariableFeatures(experiment, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(experiment)
experiment <- ScaleData(experiment, features = all.genes)
experiment <- RunPCA(experiment, features = VariableFeatures(object = experiment))
experiment <- FindNeighbors(experiment, dims = 1:10)
experiment <- FindClusters(experiment, resolution = 0.5)
experiment <- RunUMAP(experiment, dims = 1:10)

# Make the cells to start with the 
View(as.data.frame(Cells(experiment)))
View(as.data.frame(Cells(da094)))
View(as.data.frame(Cells(seq098_2)))

write.csv(Embeddings(experiment, reduction = "umap"), file = paste(outpath, '/', experiment_name, "_cell_embeddings.csv", sep=''))

colours <- colorRampPalette(brewer.pal(8, "Accent"))(length(unique(experiment$seurat_clusters)))
names(colours) <- as.character(unique(experiment$seurat_clusters))
experiment@meta.data$cellIds <- rownames(experiment@meta.data)
experiment_colours <- merge(data.frame(seurat_clusters=unique(experiment$seurat_clusters)), data.frame(cluster_colours = colours, seurat_clusters = names(colours)), by='seurat_clusters')
experiment@meta.data <- experiment@meta.data[,which(!startsWith(colnames(experiment@meta.data), 'cluster_colours'))]
experiment@meta.data <- merge(experiment@meta.data, experiment_colours, by='seurat_clusters')
rownames(experiment@meta.data) <- experiment@meta.data$cellIds
write.csv(experiment@meta.data[,c('seurat_clusters', 'cluster_colours'), drop=F], file = paste(outpath, '/', experiment_name, "_clusters.csv", sep=''))

df <- as.data.frame(experiment$seurat_clusters)
colnames(df) <- 'clusters'


for(i in 1:length(ids)){
  sampleid <- ids[i]
  print(sampleid)
  print(table(experiment$orig.ident))
  sample <- subset(experiment, subset = orig.ident == sampleid)

  df <- as.data.frame(sample$seurat_clusters)
  colnames(df) <- 'clusters'

  for (k in 1:length(unique(df$clusters))){
    cluster <- unique(df$clusters)[k]
    df.cluster <- subset(df, df$clusters == cluster)
    df.cluster$barcode <- rownames(df.cluster)

    dir.create(outpath)
    file_name <- paste(outpath, '/', sampleid, '_merged.clusters_', cluster, '.tsv', sep = '')
    print(file_name)
    write.table(sapply(str_split(df.cluster$barcode, "_"), `[[`, 1), file = file_name, row.names = F, col.names = F, quote = F)
  }

}

save(experiment, file=paste(outpath, experiment_name, ".RData", sep=''))
