# trusTEr - A trusting TE analysis
WARNING : If you are seeing this, this is a website under development. You probably are in the wrong place!

Version 0.1

Takes fastq files from 10x single cell RNA sequencing, clusters cells using Seurat, and can be used to produce 
read count matrices in a cluster level. You can also quantify reads per cluster having predefined clusters.


## Requirements
TrusTEr depends on several external software. We provide a Docker container and a conda environment for a quick-start. 

TrusTEr requires:

* Cellranger
* R (version 3.6)
    * Seurat
* TEtranscripts
* STAR aligner
* subset-bam and bamtofastq from 10x Genomics
* Velocyto

The package has been tested in Unix systems only and supports only SLURM job submissions.


## How to install 
#### Just the modules
If you fulfill the requirements, you can simply install via pip:

`pip install truster`

#### With Docker container

#### With conda environment


## Structure

trusTEr uses composition assiciation to relate three main classes: 

* Experiment: Includes information about the experiment as a whole. This is the main object you will be working with.
    * Name is required, description is optional. 
    * Register samples by providing a path to its fastq files.
* Sample: Created by giving a path to fastq files
    * Name and ID required. 
* Cluster: Created by running `getClusters` or `mergeSamples` functions.
