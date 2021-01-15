
# A workflow example

Let's begin with a simple example. Let's say you have two samples of the same tissue, they were sequenced in the same sequencing run and you want to get transposon expression per cell type found in the tissue.

<br/>
## Registering samples
The first step is to have the fastq files in a file system such as:

```
data
│
└───sample1
│   │   sample1Name_L001_I1_001.fastq.gz
│   │   sample1Name_L001_R1_001.fastq.gz
│   │   sample1Name_L001_R2_001.fastq.gz 
│   │   sample1Name_L002_I1_001.fastq.gz 
│   │   sample1Name_L002_R1_001.fastq.gz 
│   └─  sample1Name_L002_R2_001.fastq.gz 
└───sample2
    │   sample2Name_L001_I1_001.fastq.gz
    │   sample2Name_L001_R1_001.fastq.gz
    │   sample2Name_L001_R2_001.fastq.gz 
    │   sample2Name_L002_I1_001.fastq.gz 
    │   sample2Name_L002_R1_001.fastq.gz 
    └─  sample2Name_L002_R2_001.fastq.gz 
```

We can then register the samples in an experiment object

```
import truster
raw_path = ["/path/to/data/"]
example = truster.experiment(name = "example")
example.registerSamplesFromPath(raw_path)
```

This will look for this structure of files in the path and create an object "sample" per subdirectory in the path, the object will have the name of the subdirectories as ID (i.e. sample1) and the name of the files (i.e. sample1Name) as sample name. The "sample" objects are contained to the object "experiment" which here we called **example**.

If there was any sample in the `data/` folder you didn't want to include, you can unregister it with 

```
example.unregisterSample(sampleId = "sample1")
```

Similarly, you can register individual samples as 

```
example.registerSample(sampleId = "sample3", sampleName = "sample3Name", rawPath = "/path/to/sample3/fastqfiles/")
```

<br/>
## Mapping and quantification
We now need to run cellranger. TrusTEr has a wrapper for it and can be used to run it for all samples registered in an experiment. In this example:

```
crIndex = '/path/to/cellranger/index'
outdir = '/output/path'
example.quantify(crIndex, outdir)
```

If you already have the output from cellranger, feel free to set the quantification outdir with

```
example.setQuantificationOutdir(outdir)
```
Note: Before continuing, don't forget to check the quality of your samples!

<br/>
## Clustering
Clustering samples is not a trivial step. It requires knowledge of the tissue being analyzed. We wrapped a very basic clustering with Seurat as a reference to make the process smoother for the user, but be aware that there are plenty of parameters that you could  tweak to get the perfect clustering for your particular data. 
 
#### Per sample
We can get a clustering of each sample by typing

```
clusters_dir = 'wherever/you/want/the/output'
example.getClustersAllSamples(clusters_dir)
```

This will create the `clusters_dir` and a subfolder per sample (named with the `sampleIds`). Each of these subfoldes will contain a tsv file per cluster found on the sample. The tsv files contain the cells barcodes that form that cluster.

You will also get an RData file in the output directory with the Seurat object of each of your samples.

If you already have a clustering of your preference, please produce the required tsv files and set the clusters directory as:

```
example.setClustersOutdir(processClustersOutdir = outdir)
```

We could add a function in the near future that takes RData with a Seurat object and produces the tsv files in the file structure we need it. 

#### Merged samples
In some experiments, such as this one where the samples are from the same tissue, it's interesting to combine the samples and get a merged clustering. This can be achieved as

```
mergedsamples_dir = 'wherever/you/want/the/output'
example.mergeSamples(mergedsamples_dir)
```

Similarly to the clusters per sample, if you already have a clustering you want to use, you could set the directory where you contain the tsv files with the cell barcodes as 

```
example.setMergeSamplesOutdir(outdir)
```

If you are going to do this for the merged clustering, we ask you to name your tsv files as `[sampleId]_merged.clusters_[cluster number].tsv`

<br/>
## Transposon expression
Once we have the cellranger output, and the clusters that we want to check for transposon expression, we can go ahead and run trusTEr as a pipeline or step by step. 

This part is a set of seven steps:

1. Extract cell barcodes from the bam files
2. Filter for unique UMIs
3. Convert bam to fastq files
4. Concatenate lanes from the step #3 output
5. Map fastq files
6. TE quantification
7. TE quantification normalization

#### As a pipeline

You can run everything at once.

You will need to provide a gene and a TE GTF (as required by TEtranscripts output. You can check their downloads here), the path to the STAR index you want to use, the mode meaning the type of clustering you are using ("merged" or "perSample"), and the output directory. 

```
output = 'wherever/you/want/the/output'
gene_gtf = 'path/to/gene.gtf'
te_gtf = 'path/to/te.gtf'
star_index = 'path/to/star/index'

example.processClusters(mode = "merged", outdir = output, geneGTF = gene_gtf, teGTF = te_gtf, starIndex = star_index)
```

#### Step by step

If you want to wait a bit and check the output of each step of the pipeline, you an of course run it step by step.

##### 1) Extract clusters

Extract cell barcodes from the bam files typing
```
example.tsvToBamClusters(mode = "merged", outdir)
```
This will create a directory inside the outdir named `tsvToBam/`. This directory will contain subdirectories one per sample, containing the bam files for each of their clusters.

##### 2) Filter UMIs

Because PCR duplication is a needed step for 10x RNA sequencing, we need to make sure there are no duplicated molecules when quantifying repetitive elements. 

We ensure the molecules are unique by keeping only reads with a unique combination of cell barcode, UMI and sequence. 

To filter duplicates in our merged clustering bam files:

```
example.filterUMIsClusters("merged", outdir)
```

This will create a directory inside the outdir named `filterUMIs/`. Similarly to `tsvToBam/` it contains a subdirectory per sample, each containing the filtered bam file.

##### 3) Convert bam to fastq 

This is just a file conversion step dependent on bamtofastq from 10x Genomics.

```
example.bamToFastqClusters("merged", outdir) 
```

Will create a directory inside the outdir named `bamToFastq/`. Each sample subdirectory has subdirectories for each cluster which contain the fastq files of the clusters in different lanes (L00[1-9]).

##### 4) Concatenate lanes

We will take the sequence fastq files (`*_R2_001.fastq.gz`) and concatenate them. This will produce a bulk file for each cluster. 

```
example.concatenateLanesClusters("merged", outdir)
```

Again, this will create a directory inside the outdir named `concatenateLanes/` which will include sample subdirectories containing the concatenated fastq files per cluster.

##### 5) Map fastq files

We will now map these files using STAR. Again, this is just a basic wrapper using TEtranscript's authors recommendations, but if you need to tweak more parameters, feel free to do this step yourself and then setting the output directory to continue.

You will need a gene GTF and a STAR index.

```
output = 'wherever/you/want/the/output'
gene_gtf = 'path/to/gene.gtf'
star_index = 'path/to/star/index'

example.mapClusters("merged", outdir = output, geneGTF = gene_gtf, starIndex = star_index)
```

##### 6) TE quantification

This step runs TEcounts in multi mode (for details, see TEtranscripts documentation). You can download the appropiate TE GTF file here.

```
output = 'wherever/you/want/the/output'
gene_gtf = 'path/to/gene.gtf'
te_gtf = 'path/to/te.gtf'

example.TEcountsClusters("merged", outdir = output, geneGTF = gene_gtf, teGTF = te_gtf)
```

The output directory now contains a subdirectory called `TEcounts/` with samples' subdirectories and each of their clusters TE counts.

##### 7) TE counts normalization

Before continuing with the downstream analysis, we need to normalize for samples' sequencing depth and cluster size. We can do this for all samples using

```
example.normalizeTECounts("merged")
```

And that's it! You can now see your final count matrix at the output directory and if you want you can use our plot_TEexpression.R script to plot a UMAP with TE subfamilies expression. 

--------------------
Note that if you need to cancel the execution of a step (or the pipeline if you decided to go for that) you will have to wait for all samples to finish the step they are at.


