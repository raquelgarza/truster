#!/bin/env Rscript

library(ggplot2)
library(optparse)
library(Seurat)
library(stringr)
library(ggpubr)
library(data.table)
library(dplyr)
library(patchwork)
set.seed(10)

# Name of a TE subfamily
# Path to RData with Seurat object to use
# Path of TEcounts melted csv file (output from normalize_TEexpression.R)
# Mode. Merged samples or individual? (merged/individual)

# rdatas <- c("/Volumes/LaCie/FetalCortex/3_mergeSamples/fetalcortex.RData","/Volumes/LaCie/Gliomas/01.01.21/3_mergeSamples/gliomas.RData", "/Volumes/LaCie/TBI/3_mergeSamples/tbi.RData")
# tes_ids <- c("L1HS:L1:LINE","L1PA2:L1:LINE","L1PA3:L1:LINE","L1PA4:L1:LINE")
# inputs <- c("/Volumes/LaCie/FetalCortex/3_mergeSamples/clusterPipeline/TEcountsNormalized/TE_normalizedValues_aggregatedByClusters_melted.csv", "/Volumes/LaCie/Gliomas/01.01.21/3_mergeSamples/clusterPipeline/TEcountsNormalized/TE_normalizedValues_aggregatedByClusters_melted.csv", "/Volumes/LaCie/TBI/3_mergeSamples/clusterPipeline/TEcountsNormalized/TE_normalizedValues_aggregatedByClusters_melted.csv")
# outdir <- "/Volumes/LaCie/FetalCortex/3_mergeSamples/clusterPipeline/TEcountsNormalized/"
# names <- c("FetalCortex","Glioma","TBI")
# modes <- c("merged","merged","merged")
# plot_TEexpression.R -r ../3_mergedSamples/gliomas.RData -m merged -n Gliomas -t L1HS:L1:LINE,L1PA2:L1:LINE,L1PA3:L1:LINE,L1PA4:L1:LINE,L1PA5:L1:LINE,L1PA6:L1:LINE,L1PA7:L1:LINE,L1PA8:L1:LINE -i /projects/fs5/raquelgg/Gliomas/Seq073_Seq091/3_mergedSamples/clusterPipeline/TEcountsNormalized -o /projects/fs5/raquelgg/Gliomas/Seq073_Seq091/3_mergedSamples/clusterPipeline/TEplots
option_list = list(
  make_option(c("-r", "--RDatas"), type="character", default=NULL,
              help="Path to RDatas with Seurat objects (comma-delimited list)", metavar="character"),
  make_option(c("-m", "--modes"), type="character", default=NULL,
              help="Are these merged seurat objects? or individual samples? (merged/individual). In the same order of the --RDatas paths.", metavar="character"),
  make_option(c("-n", "--nameRDatas"), type="character", default=NULL,
              help="Name for experiments/samples Seurat objects (comma-delimited list). In the same order of the --RDatas paths.", metavar="character"),
  make_option(c("-t", "--teSubfamily"), type="character", default=NULL,
              help="Name of a TE subfamily or comma-delimited list", metavar="character"),
  make_option(c("-i", "--inputs"), type="character", default=NULL,
              help="Path of TEcounts melted csv files (*melted.csv output from normalize_TEexpression.R). In the same order of the --RDatas paths.", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Output directory where to place the plots produced (.pdf)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$modes) | is.null(opt$RDatas) | is.null(opt$nameRDatas) | is.null(opt$teSubfamily) | is.null(opt$inputs) | is.null(opt$outdir)){
  print_help(opt_parser)
  stop("All argument must be supplied.", call.=FALSE)
}

modes <- trimws(unlist(str_split(opt$modes, ',')))
rdatas <- trimws(unlist(str_split(opt$RDatas, ',')))
names <- trimws(unlist(str_split(opt$nameRDatas, ',')))
tes_ids <- trimws(unlist(str_split(opt$teSubfamily, ',')))
inputs <- trimws(unlist(str_split(opt$inputs, ',')))
outdir <- ifelse(endsWith(opt$outdir, "/"), opt$outdir, paste(opt$outdir, '/', sep=''))


print(c("Mode: ", modes))
print(c("RData to be used: ", rdatas))
print(c("Names to be used: ", names))
print(c("Input file: ", inputs))
print(c("Output path: ", outdir))

map2color <- function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

mypal <- colorRampPalette( c( "lightgrey", "red" ) )( 10 )

counts <- list()
for(i in 1:length(inputs)){
  input <- inputs[i]
  name <- names[i]
  counts[[name]] <- fread(input, data.table = F)  
  counts[[name]]$condition <- name
}

seurat.objs <- list()
for(i in 1:length(rdatas)){
  rdata <- rdatas[i]
  name <- names[i]
  mode <- modes[i]
  load(rdata)
  
  if(mode == 'merged'){
    seurat.objs[[name]] <- experiment
  }else{
    seurat.objs[[name]] <- sample
  }
}

colourLimits_outdir <- paste("colourLimits_", paste(sapply(str_split(tes_ids, ":"), `[[`, 1), collapse = '.'), sep='')

outdir <- paste(outdir, colourLimits_outdir, '/', sep='')

if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)  
}

te <- list()
for(i in 1:length(names)){
  name <- names[i]
  count <- counts[[name]]
  if(length(tes_ids) > 1){
    te[[name]] <- count[which(count$te_id %in% tes_ids),]
  }else{
    te[[name]] <- count[which(count$te_id == tes_ids),]
  }
}

te <- do.call(rbind, te)

# print(min(te$value))
# print(max(te$value))
te$colour <- map2color(te$value, pal=mypal)

plots <- list()
for(i in 1:length(tes_ids)){
  te_id <- tes_ids[i]
  te_subfam <- sapply(str_split(te_id, ":"), `[[`, 1)
  print(paste("TE subfamily ", te_subfam))
  te_i <- te[which(te$te_id == te_id),]
  print("################################")
  print(te_i)
  print("################################")
  for(j in 1:length(names(seurat.objs))){
    name <- names(seurat.objs)[j]
    seurat.obj <- seurat.objs[[name]]
    # print(name)
    # print(te_id)
    te_i_name <- te_i[which(te_i$condition == name & te_i$te_id == te_id),]
    te_i_name <- te_i_name[order(as.numeric(as.character(te_i_name$cluster))),]
    # print(paste("clusters ", unique(te_i_name$condition)))
    # print(unique(te_i_name$cluster))
    # print(paste("length of colours ", length(te_i_name$colour)))
    # print(te_i_name$colour)
    plots[[paste(name, te_subfam, sep = '_')]] <- DimPlot(seurat.obj, reduction='umap', cols = te_i_name$colour) 
  }
  
  title <- paste(te_subfam, paste(names(seurat.objs), collapse = '.'), sep = "_")
  
  num_of_plots <- length(paste(names(seurat.objs), te_subfam, sep = '_'))
  
  if(is.integer(num_of_plots/2)){
    wid <- 14
    hei <- 14
  }else{
    wid <- 14
    hei <- 7
  }
  
  pdf(paste(outdir, title, "_umap.pdf", sep=''), onefile = T, width = wid, height = hei)
  for (j in length(plots)) {
    # print(plots[[names(plots)[j]]])  
    print(ggarrange(plotlist=plots[paste(names(seurat.objs), te_subfam, sep = '_')],
                    labels = paste(names(seurat.objs), te_subfam, "expression"), ncol = 2))
  }
  dev.off()
  
}
# 
# te$group <- paste(te$condition, te$cluster, sep='_')
# te_unmelt <- reshape2::dcast(data = te, formula = te_id~group, value.var = "value")
# rownames(te_unmelt) <- te_unmelt$te_id
# te_unmelt <- te_unmelt[,-1]
# 
# library(pheatmap)
# pheatmap(log2(te_unmelt+0.5), cluster_cols = F)
# 


# Rscript r_scripts/plot_TEexpression.R -r /Volumes/LaCie/FetalCortex/3_mergeSamples/fetalcortex.RData,/Volumes/LaCie/Gliomas/01.01.21/3_mergeSamples/gliomas.RData,/Volumes/LaCie/TBI/3_mergeSamples/tbi.RData -m merged,merged,merged -n fetalcortex,gliomas,tbi -t L1HS:L1:LINE,L1PA2:L1:LINE,L1PA3:L1:LINE,L1PA4:L1:LINE,L1PA5:L1:LINE,L1PA6:L1:LINE -i /Volumes/LaCie/FetalCortex/3_mergeSamples/clusterPipeline/TEcountsNormalized/TE_normalizedValues_aggregatedByClusters_melted.csv,/Volumes/LaCie/Gliomas/01.01.21/3_mergeSamples/clusterPipeline/TEcountsNormalized/TE_normalizedValues_aggregatedByClusters_melted.csv,/Volumes/LaCie/TBI/3_mergeSamples/clusterPipeline/TEcountsNormalized/TE_normalizedValues_aggregatedByClusters_melted.csv -o /Volumes/LaCie/TBI/3_mergeSamples/clusterPipeline/TEcountsNormalized
