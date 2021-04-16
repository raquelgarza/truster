# trusTEr - A trusting TE analysis
WARNING : This is a website under development. You probably are in the wrong place!

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


## Introduction



## Structure

trusTEr uses composition assiciation to relate three main classes: 

* Experiment: Includes information about the experiment as a whole. This is the main object you will be working with.
    * Name is required, description is optional. 
    * Register samples by providing a path to its fastq files.
* Sample: Created by giving a path to fastq files
    * Name and ID required. 
* Cluster: Created by running `getClusters()` or `mergeSamples()` functions (Or `setClustersOutdir()` or `setMergeSamplesOutdir()`).

<img src="https://raw.githubusercontent.com/ra7555ga-s/truster/main/img/compositionAssociation.png" width="1000">

## Functionality and workflow

This package is meant to be run with the following workflow:

![alt text](https://github.com/ra7555ga-s/truster/blob/main/img/workflow.png?raw=True?=50)

Depending on the object type, you have access to different functions to go through these steps. 

`Experiment` is the main object one would work with. Here you will `registerSample()` or `registerSamplesFromPath()`. 

An object of type `Sample` has access to step 7 and some handy wrappers to use `cellRanger`, perform clustering with Seurat and run and plot RNA velocity. These functions can be called for all registered samples from your object of type `Experiment` (See `quantify()`, `getClustersAllSamples()`, `velocityAllSamples()`, `plotVelocityAllSamples()`).

The need for the class `Experiment` is clearer once the user wants to merge samples (See `mergeSamples()`) or to run the same workflow for all the samples' clusters.

The user won't work directly with an object of type `Cluster`, but this class includes all the functions needed to go through steps 1-6 of the workflow. Instead of running this pipeline individually for each cluster, one can run the workflow for each cluster of each registered sample using the `Experiment` function `processClusters(mode = "perSample", ...)`.

One can also partition the workflow and run step by step in all registered samples or in a combined clustering using the transitioning functions of `Experiment` to call the needed functions in the class `Cluster` (See `tsvToBamClusters()`, `filterUMIsClusters()`, `bamToFastqClusters()`, `concatenateLanesClusters()`, `mergeClusters()`, `mapClusters()`, `TEcountsClusters()`, `normalizeTECounts()`).


