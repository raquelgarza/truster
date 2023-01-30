#!/bin/env Rscript

library(optparse)
library(Seurat)
library(stringr)
library(data.table)
library(dplyr)
library(patchwork)
set.seed(10)
# 
# mode <- "merged" 
# group_name <- "control"
# by_factor <- "seurat_clusters" 
# samples <- c("DA428","DA488","ASAP14","DA480","DA806","DA795","DA785","DA804_ASAP49_Ctl_NP18-159_PFCTX_Seq163A")
# outdir <- "/Volumes/My Passport/TBI/03.05.22/3_combinedUMAP_perCluster_ASAP_controls/clusterPipeline/TE_counts_normalized/multiple/"
# indir <- "/Volumes/My Passport/TBI/03.05.22/3_combinedUMAP_perCluster_ASAP_controls/clusterPipeline/TE_counts/multiple/"
# rds <- "/Volumes/My Passport/TBI/03.05.22/3_combinedUMAP_perCluster_ASAP_controls/tbi_ASAP_controls.rds"
# obj_name <- "tbi_ASAP_controls"

# Bed file of transposons
# Tab file of classification of transposons
# Path to RData 
# Path to TEcounts folder with output per cluster

option_list = list(
  make_option(c("-m", "--mode"), type="character", default=NULL,
              help="Merged samples or individual? (merged/individual)", metavar="character"),
  make_option(c("-r", "--RDS"), type="character", default=NULL,
              help="Path to RDS with Seurat object", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="Sample id or name", metavar="character"),
  make_option(c("-g", "--groupName"), type="character", default=NULL,
              help="Group name", metavar="character"),
  make_option(c("-f", "--byFactor"), type="character", default="seurat_clusters",
              help="Groupping factor (Default = seurat_clusters)", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="Sample ids in the group (comma delimited)", metavar="character"),
  make_option(c("-i", "--indir"), type="character", default=NULL,
              help="Path to TEcounts folder with output per cluster or parent directory from all samples folder", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Output directory for the matrix with normalized counts", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$indir) | is.null(opt$outdir) | is.null(opt$RDS) | is.null(opt$mode) ){
  print_help(opt_parser)
  stop("Indir, outdir, rds and mode arguments must be supplied.", call.=FALSE)
}

mode <- trimws(opt$mode)
rds <- opt$RDS
obj_name <- opt$name
group_name <- opt$groupName
by_factor <- if(grepl(",", opt$byFactor)) unlist(str_split(opt$byFactor, ",")) else opt$byFactor
samples <- unlist(str_split(opt$samples, ","))
indir <- ifelse(endsWith(opt$indir, "/"), opt$indir, paste(opt$indir, '/', sep=''))
outdir <- ifelse(endsWith(opt$outdir, "/"), opt$outdir, paste(opt$outdir, '/', sep=''))

if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

print(c("Mode: ", mode))
print(c("RDS to be used: ", rds))
print(c("Input path: ", indir))
print(c("Output path: ", outdir))

# Load object
seurat.obj <- readRDS(rds)
# Keeping only the samples that are supposed to be included (only in case of it being more than one sample in the object)
if(mode == "merged"){
  seurat.obj <- subset(seurat.obj, orig.ident %in% samples)  
}

# If there are more than one grouping factor, create one metadata column for the different combinations of it. If there is just one factor, the metadata column will be called the same
factors_by_cell <- FetchData(seurat.obj, vars = by_factor)
factors_by_cell$groupping_factor <- apply(factors_by_cell, 1, function(x) paste(x, collapse = "_"));
seurat.obj <- AddMetaData(seurat.obj, metadata = factors_by_cell$groupping_factor, col.name = "groupping_factor")

cluster_sizes <- data.frame(cluster.size=table(FetchData(seurat.obj, "groupping_factor")))  
colnames(cluster_sizes) <- c('cluster', 'cluster.size')

files <- list.files(indir, recursive = F)

if( mode == "merged"){
  cluster_sizes$name <- paste(obj_name, group_name, cluster_sizes$cluster, sep='_')
  regex_startswith <- paste("^", cluster_sizes$name, sep="", collapse = "|")
  files <- files[which(grepl(regex_startswith, files))]
}else{
  cluster_sizes$name <- paste(obj_name, cluster_sizes$cluster, sep='_')
}

coldata <- data.frame()
for(i in 1:length(files)){
  if(mode == 'merged'){
    cluster <- unlist(str_split(unlist(str_split(files[i], paste(obj_name, group_name, "", sep="_")))[2], "_.cntTable"))[1]
    name <- paste(obj_name, group_name, cluster, sep="_") 
    coldata <- rbind(coldata, data.frame(sample=obj_name, name = name, cluster=cluster)) 
  }else{
    file_name <- sub("_$", "", sapply(str_split(files[i], '.cntTable'), `[[`, 1))
    cluster <- unlist(str_split(file_name, '_'))[2]
    name <- paste(obj_name, cluster, sep="_")
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

coldata <- merge(coldata, cluster_sizes[,c('name', 'cluster.size')], by='name')
rownames(coldata) <- coldata$name

TEcounts <- subset(TEcounts, !startsWith(TEcounts$TE, 'ENS'))

#### TE_norm_cluster_size ####
rownames(coldata) <- coldata$name
te_counts <- TEcounts[, rownames(coldata)]
te_counts_size_norm <- TEcounts[, rownames(coldata)]
te_counts_size_norm[] <- mapply('/', te_counts_size_norm[, rownames(coldata)], coldata$cluster.size)

te_counts_size_norm$te_id <- rownames(te_counts_size_norm)
te_counts_size_norm_melt <- reshape2::melt(te_counts_size_norm, by=list(c('te_id')))
te_counts_size_norm_melt$cluster <- sapply(str_split(te_counts_size_norm_melt$variable, paste(obj_name, group_name, "", sep="_")),`[[`, 2)

cells_clusters <- FetchData(seurat.obj, vars = "groupping_factor")
colnames(cells_clusters) <- "cluster"
cells_clusters$cell_ids <- rownames(cells_clusters)

te_counts_size_norm_melt_percell <- merge(cells_clusters, te_counts_size_norm_melt, by='cluster', all = T)
te_counts_size_norm_melt_percell <- reshape2::dcast(te_counts_size_norm_melt_percell, formula = te_id~cell_ids, value=value)

rownames(te_counts_size_norm_melt_percell) <- te_counts_size_norm_melt_percell$te_id  
te_counts_size_norm_melt_percell <- te_counts_size_norm_melt_percell[,-1]

# Norm by cluster size
seurat.obj[["TE_norm_cluster_size"]] <- CreateAssayObject(counts = te_counts_size_norm_melt_percell)

#### TE_raw ####
te_counts$te_id <- rownames(te_counts)
te_counts_melt <- reshape2::melt(te_counts, by=list(c('te_id')))
te_counts_melt$cluster <- sapply(str_split(te_counts_melt$variable, paste(obj_name, group_name, "", sep="_")),`[[`, 2)
# te_counts_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_melt$variable, paste("_", group_name, "_", sep="")),`[[`, 2))

te_counts_melt_percell <- merge(cells_clusters, te_counts_melt, by='cluster', all = T)
te_counts_melt_percell <- reshape2::dcast(te_counts_melt_percell, formula = te_id~cell_ids, value=value)
rownames(te_counts_melt_percell) <- te_counts_melt_percell$te_id  
te_counts_melt_percell <- te_counts_melt_percell[,-1]

seurat.obj[["TE_raw"]] <- CreateAssayObject(counts = te_counts_melt_percell)

#### TE_norm_cluster_size_num_reads ####
te_counts <- TEcounts[, rownames(coldata)]
te_counts_size_norm <- TEcounts[, rownames(coldata)]
te_counts_size_norm[] <- mapply('/', te_counts_size_norm[, rownames(coldata)], coldata$cluster.size)

tmp <- FetchData(seurat.obj, vars=c("groupping_factor") )
raw_counts <- GetAssayData(seurat.obj, slot = "counts", assay="RNA")
total_num_reads <- sum(colSums(raw_counts))

te_counts_size_norm$te_id <- rownames(te_counts_size_norm)
te_counts_size_reads_norm <- te_counts_size_norm[, rownames(coldata)] / total_num_reads[rownames(coldata)]
te_counts_size_reads_norm <- log1p(te_counts_size_reads_norm[, rownames(coldata)] * 1e+7)

te_counts_size_reads_norm$te_id <- rownames(te_counts_size_reads_norm)
te_counts_size_reads_norm_melt <- reshape2::melt(te_counts_size_reads_norm, by=list(c('te_id')))
te_counts_size_reads_norm_melt$cluster <- sapply(str_split(te_counts_size_reads_norm_melt$variable, paste(obj_name, group_name, "", sep="_")),`[[`, 2)
# te_counts_size_reads_norm_melt$cluster <- gsub("cluster_", "", sapply(str_split(te_counts_size_reads_norm_melt$variable, paste("_", group_name, "_", sep="")),`[[`, 2))

te_counts_size_reads_norm_melt_percell <- merge(cells_clusters, te_counts_size_reads_norm_melt, by='cluster', all = T)
te_counts_size_reads_norm_melt_percell <- reshape2::dcast(te_counts_size_reads_norm_melt_percell, formula = te_id~cell_ids, value=value)
rownames(te_counts_size_reads_norm_melt_percell) <- te_counts_size_reads_norm_melt_percell$te_id  
te_counts_size_reads_norm_melt_percell <- te_counts_size_reads_norm_melt_percell[,-1]

# Norm by cluster size and sequencing depth
seurat.obj[["TE_norm_cluster_size_num_reads"]] <- CreateAssayObject(counts = te_counts_size_reads_norm_melt_percell)

te_counts_size_norm <- te_counts_size_norm[,c(colnames(te_counts_size_norm)[ncol(te_counts_size_norm)], colnames(te_counts_size_norm)[-ncol(te_counts_size_norm)])]
te_counts_size_reads_norm <- te_counts_size_reads_norm[,c(colnames(te_counts_size_reads_norm)[ncol(te_counts_size_reads_norm)], colnames(te_counts_size_reads_norm)[-ncol(te_counts_size_reads_norm)])]
te_counts <- te_counts[,c(colnames(te_counts)[ncol(te_counts)], colnames(te_counts)[-ncol(te_counts)])]
te_counts$te_id <- rownames(te_counts)
te_counts <- te_counts[,c("te_id",colnames(te_counts)[-ncol(te_counts)])]

file_name <- ifelse(is.null(group_name), "TE_norm_cluster_size_matrix.csv", paste(group_name, "TE_norm_cluster_size_matrix.csv", sep="_"))
fwrite(te_counts_size_norm, paste(outdir, file_name, sep=''), quote = F, row.names = F, col.names = T, verbose = T)
file_name <- ifelse(is.null(group_name), "TE_norm_cluster_size_num_reads_matrix.csv", paste(group_name, "TE_norm_cluster_size_num_reads_matrix.csv", sep="_"))
fwrite(te_counts_size_reads_norm, paste(outdir, file_name, sep=''), quote = F, row.names = F, col.names = T, verbose = T)
# To compare with bulk data
file_name <- ifelse(is.null(group_name), "TE_raw_matrix.csv", paste(group_name, "TE_raw_matrix.csv", sep="_"))
fwrite(te_counts, paste(outdir, file_name, sep=''), quote = F, row.names = F, col.names = T, verbose = T)

rds <- paste(unlist(str_split(rds, ".rds"))[1], "_", group_name, ".rds", sep="")
saveRDS(seurat.obj, file = rds)



