#!/usr/bin/python3

import anndata
import scvelo as scv
import pandas as pd
import re
import os
import subprocess
import sys
import getopt

def main(argv):
    loom_file = ''
    umap_file = ''
    clusters_file = ''
    sample_name = ''
    outdir = ''
    try:
        opts, args = getopt.getopt(argv,"hl:n:u:c:",["loom=", "sample_name=", "umap=", "clusters=", "outdir="])
    except getopt.GetoptError:
        print('plotVelocity -l <loom> -n <sample_name> -u <umap> -c <clusters> -o <outdir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('plotVelocity -l <loom> -n <sample_name> -u <umap> -c <clusters> -o <outdir>')
            sys.exit()
        elif opt in ("-l", "--loom"):
            loom_file = arg
        elif opt in ("-n", "--sample_name"):
            sample_name = arg
        elif opt in ("-u", "--umap"):
            umap_file = arg
        elif opt in ("-c", "--clusters"):
            clusters_file = arg
        elif opt in ("-o", "--outdir"):
            outdir = arg

    # >>> sample = anndata.read_loom("FetalCortex/recovered_Seq030.Seq041.Seq044_FetalCortex/12.10.20_premRNA/1_counts/DA094/velocyto/DA094.loom")
    # Variable names are not unique. To make them unique, call `.var_names_make_unique`.
    # >>> cell_clusters = pd.read_csv("FetalCortex/recovered_Seq030.Seq041.Seq044_FetalCortex/12.10.20_premRNA/5_rnavelocity/DA094_clusters.csv")
    # >>> sample_umap = pd.read_csv("FetalCortex/recovered_Seq030.Seq041.Seq044_FetalCortex/12.10.20_premRNA/5_rnavelocity/DA094_cell_embeddings.csv")
    sample = anndata.read_loom(loom_file)
    sample_umap = pd.read_csv(umap_file)
    cell_clusters = pd.read_csv(clusters_file)

    sample_index = pd.DataFrame(sample.obs.index)
    sample_index = sample_index.rename(columns = {0:'Cell ID'})

    sample_umap = sample_umap.rename(columns = {'Unnamed: 0':'Cell ID'})
    cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})

    cell_clusters['Cell ID'] = [
            '%s:%sx' % (cluster, re.sub('-[0-9]+$', '', cell))
            for cluster, cell in zip([sample_name] * cell_clusters.shape[0], cell_clusters['Cell ID'].tolist())
        ]

    sample_umap['Cell ID'] = [
            '%s:%sx' % (cluster, re.sub('-[0-9]+$', '', cell))
            for cluster, cell in zip([sample_name] * sample_umap.shape[0], sample_umap['Cell ID'].tolist())
        ]
    keep_cell = [cell in sample_umap['Cell ID'].tolist() for cell in sample_index['Cell ID'].tolist()]
    sample = sample[keep_cell]

    umap_ordered = sample_index.merge(sample_umap,on="Cell ID")
    umap_ordered = umap_ordered.iloc[:,1:]
    sample.obsm['X_umap'] = umap_ordered.values

    sample.obs['Cell ID'] = sample.obs.index
    sample.obs = sample.obs.merge(cell_clusters,on="Cell ID")
    sample.obs.index = sample.obs['Cell ID']

    scv.pp.filter_and_normalize(sample)
    scv.pp.moments(sample)
    scv.tl.velocity(sample, mode = "stochastic")
    scv.tl.velocity_graph(sample)

    colours = cell_clusters[["seurat_clusters", "cluster_colours"]].drop_duplicates()
    colours = {str(cluster):colour for cluster,colour in zip(colours.seurat_clusters, colours.cluster_colours)}

    scv.pl.velocity_embedding(sample, arrow_length=5, arrow_size=3, dpi=120, color='seurat_clusters', palette=colours, save = (sample_name + "_embedding.pdf"), show=False)
    subprocess.call(["mv", ("./figures/scvelo_" + sample_name + "_embedding.pdf"), outdir])

    scv.pl.velocity_embedding_stream(sample, basis='umap', color='seurat_clusters', palette=colours, save = (sample_name + "_stream.pdf"), show=False)
    subprocess.call(["mv", ("./figures/scvelo_" + sample_name + "_stream.pdf"), outdir])

    scv.pl.umap(sample, color='seurat_clusters', save = (sample_name + "_UMAP.pdf"), show=False)
    subprocess.call(["mv", ("./figures/scvelo_" + sample_name + "_UMAP.pdf"), outdir])
    #scv.pl.proportions(sample, groupby='seurat_clusters')

if __name__ == "__main__":
   main(sys.argv[1:])