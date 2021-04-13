import subprocess
import json
from .jobHandler import *
import os
from .sample import Sample
from .cluster import Cluster
from .bcolors import Bcolors
import concurrent.futures
import copy

class Experiment:

    def __init__(self, name="", slurmPath=None, modulesPath=None):
        self.name = name
        self.slurmPath = slurmPath
        self.modulesPath = modulesPath
        self.samples = {}
        self.logfile = self.name + ".log"

        with open(self.logfile, "w+") as log:
            msg = "Project " + self.name + " created.\n"
            log.write(msg)
            if slurmPath != None:
                try:
                    with open(slurmPath, "r") as config_file:
                        self.slurm = json.load(config_file)
                        msg = "Configuration file loaded.\n"
                        log.write(msg)
                except FileNotFoundError:
                    msg = Bcolors.FAIL + "Error: Cluster configuration file not found" + Bcolors.ENDC + "\n"
                    print(msg)
                    log.write(msg)
                    return 1

            if modulesPath != None:
                try:
                    with open(modulesPath, "r") as modules_file:
                        self.modules = json.load(modules_file)
                        msg = "Software modules json loaded.\n"
                        log.write(msg)
                except FileNotFoundError:
                    msg = Bcolors.FAIL + "Error: Module configuration file not found" + Bcolors.ENDC + "\n"
                    print(msg)
                    log.write(msg)
                    return 1

    def registerSample(self, sampleId = "", sampleName = "", rawPath = ""):
        newSample = {sampleId : Sample(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = rawPath, logfile = self.logfile)}
        self.samples = {**self.samples, **newSample}
        
        with open(self.logfile, "a") as log:
            msg = "Sample " + sampleId + " registered.\n"
            log.write(msg)
        # self.samples.extend(sampleTruster(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName))

    def unregisterSample(self, sampleId):
        self.samples.pop(sampleId)
        with open(self.logfile, "a") as log:
            msg = "Sample " + sampleId + " unregistered.\n"
            log.write(msg)

    # https://florimond.dev/blog/articles/2018/08/python-mutable-defaults-are-the-source-of-all-evil/
    # Change indir = []
    def registerSamplesFromPath(self, indir=[], folderNamesAsSampleIds=True):
        if(folderNamesAsSampleIds):
            for i in indir:
                samples = {item[0] : item[1].split("_L00")[0] for item in [(os.path.join(dp, f)).split("/")[-2:] for dp, dn, fn in os.walk(os.path.expanduser(i)) for f in fn]}
                newSamples = {sampleId : Sample(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = i, logfile = self.logfile) for sampleId, sampleName in samples.items()}
                self.samples = {**self.samples, **newSamples}
                # self.samples.extend([sampleTruster(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = i) for sampleId, sampleName in samples.items()])
            with open(self.logfile, "a") as log:
                msg = "Registered samples: " + str(', '.join([sample.sampleId for sample in list(self.samples.values())]) + ".\n")
                log.write(msg)

    def quantify(self, crIndex, outdir, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:

                for sample in self.samples.values():
                    sampleIndir = os.path.join(sample.rawPath, sample.sampleId)
                    sampleOutdir = os.path.join(outdir, sample.sampleId)
                    executor.submit(sample.quantify, crIndex, sampleIndir, sampleOutdir)
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + ".\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    def setQuantificationOutdir(self, sampleId, cellranger_outdir):
        self.samples[sampleId].setQuantificationOutdir(cellranger_outdir)

    def getClustersAllSamples(self, outdir, res = 0.5, percMitochondrial = None, minGenes = None, maxGenes = None, normalizationMethod = "LogNormalize", excludeFilesPath=None, maxSize = 500, jobs=1):
        with open(self.logfile, "a") as log:
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample in list(self.samples.values()):
                        if os.path.isdir(sample.quantifyOutdir):
                            sampleIndir = sample.quantifyOutdir
                        else:
                            msg = "Error: File not found. Please make sure that " + sample.quantifyOutdir + " exists.\n"
                            log.write(msg)
                            return 1
                        sampleOutdir = os.path.join(outdir, sample.sampleId)
                        res = str(res)
                        maxSize = str(maxSize)
                        minGenes = str(minGenes)
                        maxGenes = str(maxGenes)

                        if excludeFilesPath != None:
                            if os.path.isfile(os.path.join(excludeFilesPath, (sample.sampleId + "_exclude.tsv"))):
                                excludeFile = os.path.join(excludeFilesPath, (sample.sampleId + "_exclude.tsv"))
                            elif os.path.isfile(os.path.join(excludeFilesPath, (sample.sampleName + "_exclude.tsv"))):
                                excludeFile = os.path.join(excludeFilesPath, (sample.sampleName + "_exclude.tsv"))
                            else:
                                msg = "Error: File not found. Please make sure that " + excludeFilesPath + " contains a file named as the sample ID or name. E.g. " + os.path.join(excludeFilesPath, (sample.sampleName + "_exclude.tsv")) + " or " + os.path.join(excludeFilesPath, (sample.sampleId + "_exclude.tsv"))
                                log.write(msg)
                                return 1
                            
                            if os.path.isfile(excludeFile):
                                msg = "Clustering " + sample.sampleId + " excluding from " + excludeFile + "\n"
                                executor.submit(sample.getClusters, sampleIndir, sampleOutdir, res, percMitochondrial, minGenes, maxGenes, normalizationMethod, excludeFile, maxSize)
                            else:
                                msg = "Clustering " + sample.sampleId + " using all cells. File " + excludeFile + " not found.\n"
                                executor.submit(sample.getClusters, sampleIndir, sampleOutdir, res, percMitochondrial, minGenes, maxGenes, normalizationMethod, None, maxSize)
                            self.clustersExclusiveOutdir = outdir
                        else:
                            msg = "Clustering " + sample.sampleId + " using all cells.\n"
                            executor.submit(sample.getClusters, sampleIndir, sampleOutdir, res, percMitochondrial, minGenes, maxGenes, normalizationMethod, None, maxSize)
                            self.clustersOutdir = outdir
                        log.write(msg)
                
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + ".\n"
                log.write(msg)

    def setClustersOutdir(self, clustersOutdir):
        with open(self.logfile, "a") as log:
            self.clustersOutdir = clustersOutdir
            # Register clusters in each sample 
            for sample in list(self.samples.values()):
                for i in os.listdir(clustersOutdir):
                    if os.path.isdir(os.path.join(clustersOutdir, i)) and  i == sample.sampleId:
                        msg = sample.registerClustersFromPath(os.path.join(clustersOutdir, sample.sampleId)) 
                        log.write(msg)
                    elif os.path.isdir(os.path.join(clustersOutdir, i)) and  i == sample.sampleName:
                        msg = sample.registerClustersFromPath(os.path.join(clustersOutdir, sample.sampleId)) 
                        log.write(msg)
            msg = "The directory for clusters of individual samples is set to: " + clustersOutdir + ".\n"
            log.write(msg)

    def setClustersExclusiveOutdir(self, clustersExclusiveOutdir):
        with open(self.logfile, "a") as log:
            self.clustersExclusiveOutdir = clustersExclusiveOutdir
            # Register clusters in each sample 
            for sample in list(self.samples.values()):
                for i in os.listdir(clustersExclusiveOutdir):
                    if os.path.isdir(os.path.join(clustersExclusiveOutdir, i)) and  i == sample.sampleId:
                        msg = sample.registerClustersFromPath(os.path.join(clustersExclusiveOutdir, sample.sampleId))
                        log.write(msg)
                    elif os.path.isdir(os.path.join(clustersExclusiveOutdir, i)) and  i == sample.sampleName:
                        msg = sample.registerClustersFromPath(os.path.join(clustersExclusiveOutdir, sample.sampleId))
                        log.write(msg)
            msg = "The directory for clusters of individual samples is set to: " + clustersExclusiveOutdir + ". Which includes selected cells. Samples' clustering is now without the excluded cells.\n"
            log.write(msg)

    def velocityAllSamples(self, teGTF, geneGTF, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    if os.path.isdir(sample.quantifyOutdir):
                            sampleIndir = sample.quantifyOutdir
                    else:
                        msg = "Error: File not found. Please make sure that " + sample.quantifyOutdir + " exists.\n"
                        log.write(msg)
                        return 1
                    # args = [teGTF, geneGTF, sampleIndir]
                    # executor.submit(lambda p: sample.velocity(*p), args)
                    executor.submit(sample.velocity, teGTF, geneGTF, sampleIndir)
                executor.shutdown(wait=True)
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + "\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    def plotVelocityAllSamples(self, indir, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    if os.path.isdir(sample.quantifyOutdir):
                            sampleIndir = sample.quantifyOutdir
                    else:
                        msg = "Error: File not found. Please make sure that " + sample.quantifyOutdir + " exists.\n"
                        log.write(msg)
                        return 1
                    loom = os.path.join(sampleDir, "velocyto", (sample.sampleId + ".loom"))
                    sampleOutdir = os.path.join(sampleDir, "velocyto", "plots")
                    sampleIndir = os.path.join(indir, sample.sampleId)
                    executor.submit(sample.plotVelocity, loom, sampleIndir, sampleOutdir)
                executor.shutdown(wait=True)
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + "\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    def plotVelocityMerged(self, loom, names):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                cwd = os.path.dirname(os.path.realpath(__file__))
                os.path.join(cwd, "py_scripts/plotVelocity.py")
                outdir = self.outdirMergedClusters
                # print('plotVelocity -l <loom> -n <sample_name> -u <umap> -c <clusters> -o <outdir>')
                cmd = ["python", os.path.join(cwd, "py_scripts/plotVelocity"), "-l", ','.join(loom), "-n", ','.join(names), "-u", ','.join([os.path.join(outdir, (name + "_cell_embeddings.csv")) for name in names]), "-c", ','.join([os.path.join(outdir, (name + "_clusters.csv")) for name in names]), "-o", outdir]
    
                if self.slurm != None:
    
                    cmd = ' '.join(cmd)
    
                    jobFile =  os.path.join("velocity_scripts/", (self.sampleId + "_plotVelocity.sh"))
                    try:
                        jobId = runJob("plotVelocity", jobFile, cmd, self.slurm, self.modules)
                        msg = sucessSubmit("plotVelocity", self.sampleId, jobId)
                        log.write(msg)

                        exitCode = waitForJob(jobId)
                        msg = checkExitCodes("plotVelocity", ("Sample " + self.sampleId), jobId, exitCode)
                        log.write(msg)
                        
                    except:
                        msg = genericError("plotVelocity", self.sampleId)
                        log.write(msg)
                        return waitForJob(jobId)
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    #def setProcessClustersOutdir(self, processClustersOutdir):
    #    self.processedClustersOutdir = processClustersOutdir

    def mergeSamples(self, outdir, normalizationMethod, integrateSamples = "FALSE", maxSize=500):
        # Rscript {input.script} -i {rdata} -n {samplenames} -o {params.outpath}
        # Paths to RData files
        with open(self.logfile, "a") as log:
            if hasattr(self, 'clustersExclusiveOutdir'):
                samplesSeuratRdata = [os.path.join(self.clustersExclusiveOutdir, sample.sampleId, (sample.sampleId + ".rds")) for sample in list(self.samples.values())]
            else:
                samplesSeuratRdata = [os.path.join(self.clustersOutdir, sample.sampleId, (sample.sampleId + ".rds")) for sample in list(self.samples.values())]
            
            # Sample ids
            samplesIds = [sample.sampleId for sample in list(self.samples.values())]
            
            msg = "Merging samples to produce a combined clustering.\n"
            log.write(msg)
            # If we haven't made the merge before, create a directory to store the scripts needed to do so
            if not os.path.exists("mergeSamples_scripts"):
                os.makedirs("mergeSamples_scripts", exist_ok=True)
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            
            maxSize = str(maxSize)
            # Run script ../r_scripts/merge_samples.R with input (-i) of the RData paths
            # and output (-o) of the output directory desired, -s for sample ids,
            # and -e for sample names used in cellranger
            cwd = os.path.dirname(os.path.realpath(__file__))
            cmd = ["Rscript", os.path.join(cwd, "r_scripts/merge_samples.R"), "-i", ','.join(samplesSeuratRdata), "-o", outdir, "-s", ','.join(samplesIds), "-e", self.name, "-n", normalizationMethod, "-S", maxSize, "-I", integrateSamples]
            
            # If we are on a server with slurm, use the configuration file to send the jobs 
            if self.slurmPath != None:    
                cmd = ' '.join(cmd)
                
                # This command will be stored along with the slurm configurations at
                # mergeSamples_scripts/sample_mergeSamples.sh
                jobFile =  "mergeSamples_scripts/" + self.name + "_mergeSamples.sh"
                try:
                    # Run job script
                    jobId = runJob("mergeSamples", jobFile, cmd, self.slurm, self.modules)
                    
                    # Print if the job was succesfully submitted
                    msg = sucessSubmit("mergeSamples", self.name, jobId)
                    print(msg)
                    log.write(msg)
    
                    # Wait for the job to finish and returns an exit code
                    exitCode = waitForJob(jobId)
                    print("Exit code: " + str(exitCode))
                    # Print if the job was finished succesfully or not
                    msg = checkExitCodes("mergeSamples", ("Experiment " + self.name), jobId, exitCode)
                    print(msg)
                    log.write(msg)
                    
                    # If it finished succesfully then 
                    if exitCode == 0:
                
                        # For each of the registered samples
                        self.mergeSamples = copy.deepcopy(self.samples)
    
                        # Empty the clusters bc we made new ones (shared/merged)
                        for k,v in self.mergeSamples.items():
                            v.emptyClusters()
                    
                        msg = "Emptied clusters"
                        print(msg)
                        log.write(msg)
                    
                        # print([j.clusters for j in self.mergeSamples.values()])
    
                        # Make a dictionary of the same sort as the registered samples
                        # for example {sample1 : [cluster1, cluster2]}
                        # with the clusters that we created in the outdir we passed to R
                        # They all have the words "merged.clusters", after that is the number
                        # Before that is the sample id, which we can use as a key in the 
                        # mergeSamplesClusters dictionary and just append the cluster objects
                        # To the empty list we now have
                        for i in os.listdir(outdir):
                            if(i.endswith(".tsv")):
                                clusterName = i.split(".tsv")[0]
                                sampleId = clusterName.split("_merged.clusters")[0]
                                cluster = Cluster(clusterName = clusterName, tsv = os.path.join(outdir, i), logfile = self.logfile)
                                # print(sampleId, cluster)
                                # print(self.mergeSamples[sampleId].clusters)
                                self.mergeSamples[sampleId].clusters.append(cluster)
                        # print([j.clusters for j in self.mergeSamples.values()])
                        self.mergeSamplesOutdir = outdir
                    return exitCode
                except:
                    msg = genericError("mergeSamples", self.name)
                    print(msg)
                    log.write(msg)
                    return
            else:
                subprocess.call(cmd)
                self.mergeSamplesOutdir = outdir

    def setMergeSamplesOutdir(self, mergeSamplesOutdir):
        self.mergeSamplesOutdir = mergeSamplesOutdir
        # For each of the registered samples
        self.mergeSamples = copy.deepcopy(self.samples)
    
        # Empty the clusters bc we made new ones (shared/merged)
        for k,v in self.mergeSamples.items():
            v.emptyClusters()

        for i in os.listdir(mergeSamplesOutdir):
            if(i.endswith(".tsv")):
                clusterName = i.split(".tsv")[0]
                sampleId = clusterName.split("_merged.clusters")[0]
                cluster = Cluster(clusterName = clusterName, tsv = os.path.join(mergeSamplesOutdir, i), logfile = self.logfile)

                # print(sampleId, cluster)
                # print(self.mergeSamples[sampleId].clusters)
                self.mergeSamples[sampleId].clusters.append(cluster)

        with open(self.logfile, "a") as log:
            msg = "The directory for clusters of combined samples is set to: " + mergeSamplesOutdir + ".\n\n"
            log.write(msg)

            registeredClusters = [sample.clusters for sampleId,sample in self.mergeSamples.items()]
            namesRegisteredClusters = ', '.join([cluster.clusterName for listOfClusters in registeredClusters for cluster in listOfClusters])
            msg = "Registered mergeSamples clusters per sample: " + namesRegisteredClusters + "\n\n"
            log.write(msg)

    def tsvToBamClusters(self, mode, outdir, jobs=1):
        print("Running tsvToBam with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samplesDict = self.mergeSamples
                else:
                    if mode == "perSample":
                        samplesDict = self.samples
                    else:
                        msg = "Please specify a mode (merged/perSample).\n"
                        print(msg)
                        log.write(msg)
                        return 2
                    
                self.tsvToBam_results = []
                msg = "Extracting cell barcodes from BAM files.\n"
                log.write(msg)

                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            if os.path.isdir(sample.quantifyOutdir):
                                bam = os.path.join(sample.quantifyOutdir, "outs/possorted_genome_bam.bam")
                            else:
                                msg = "Error: File not found. Please make sure that " + sample.quantifyOutdir + " exists.\n"
                                log.write(msg)
                                return 1
                            outdir_sample = os.path.join(outdir, "tsvToBam/", sampleId)
                            self.tsvToBam_results.append(executor.submit(cluster.tsvToBam, sampleId, bam, outdir_sample, self.slurm, self.modules))

                tsvToBam_exitCodes = [i.result()[1] for i in self.tsvToBam_results]
                tsvToBam_allSuccess = all(exitCode == 0 for exitCode in tsvToBam_exitCodes)
                
                if tsvToBam_allSuccess:
                    msg = "\ntsvToBam finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\ntsvToBam did not finished succesfully for all samples\n"
                    log.write(msg)
                    return False

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing tsvToBam for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def filterUMIsClusters(self, mode, outdir, jobs=1):
        print("Running filterUMIs with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samplesDict = self.mergeSamples
                else:
                    if mode == "perSample":
                        samplesDict = self.samples
                    else:
                        msg = "Please specify a mode (merged/perSample).\n"
                        print(msg)
                        log.write(msg)
                        return 2

                self.filterUMIs_results = []
                msg = "Extracting cell barcodes from BAM files.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            inbam = os.path.join(outdir, "tsvToBam/", sampleId, (cluster.clusterName + ".bam"))
                            outdir_sample = os.path.join(outdir, "filterUMIs/", sampleId)
                            self.filterUMIs_results.append(executor.submit(cluster.filterUMIs, sampleId, inbam, outdir_sample, self.slurm, self.modules))
                            
                filterUMIs_exitCodes = [i.result()[1] for i in self.filterUMIs_results]
                filterUMIs_allSuccess = all(exitCode == 0 for exitCode in filterUMIs_exitCodes)
                
                if filterUMIs_allSuccess:
                    msg = "\nfilterUMIs finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nfilterUMIs did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing filterUMIs for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)
 
    def bamToFastqClusters(self, mode, outdir, jobs=1):
        print("Running bamToFastq with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samplesDict = self.mergeSamples
                else:
                    if mode == "perSample":
                        samplesDict = self.samples
                    else:
                        msg = "Please specify a mode (merged/perSample).\n"
                        print(msg)
                        log.write(msg)
                        return 2

                self.bamToFastq_results = []
                msg = "Converting BAM to FastQ files.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            bam = os.path.join(outdir, "filterUMIs/", sampleId, (cluster.clusterName + "_filtered.bam"))
                            outdir_sample = os.path.join(outdir, "bamToFastq/", sampleId)
                            self.bamToFastq_results.append(executor.submit(cluster.bamToFastq, sampleId, bam, outdir_sample, self.slurm, self.modules))
                            
                bamToFastq_exitCodes = [i.result()[1] for i in self.bamToFastq_results]
                bamToFastq_allSuccess = all(exitCode == 0 for exitCode in bamToFastq_exitCodes)

                if bamToFastq_allSuccess:
                    msg = "\nbamToFastq finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nbamToFastq did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing bamToFastq for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def concatenateLanesClusters(self, mode, outdir, jobs=1):
        print("Running concatenateLanes with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samplesDict = self.mergeSamples
                else:
                    if mode == "perSample":
                        samplesDict = self.samples
                    else:
                        msg = "Please specify a mode (merged/perSample).\n"
                        print(msg)
                        log.write(msg)
                        return

                self.concatenateLanes_results = []
                msg = "Concatenating FastQ files.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            indir = os.path.join(outdir, "bamToFastq/", sampleId)
                            outdir_sample = os.path.join(outdir, "concatenateLanes/", sampleId)
                            self.concatenateLanes_results.append(executor.submit(cluster.concatenateLanes, sampleId, indir, outdir_sample, self.slurm, self.modules))
                            
                concatenateLanes_exitCodes = [i.result()[1] for i in self.concatenateLanes_results]
                concatenateLanes_allSuccess = all(exitCode == 0 for exitCode in concatenateLanes_exitCodes)

                if concatenateLanes_allSuccess:
                    msg = "\nConcatenateLanes finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nConcatenateLanes did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing concatenateLanes for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def mergeClusters(self, outdir):
        with open(self.logfile, "a") as log:
            try:
                outdirMergedClusters = os.path.join(outdir, "mergedClusters")
                if not os.path.exists(outdirMergedClusters):
                    os.makedirs(outdirMergedClusters, exist_ok=True)

                mergeSamplesListsClusters = [samplesClusters.clusters for samplesClusters in self.mergeSamples.values()]
                mergeSamplesClusters = [cluster.clusterName.split("merged.clusters_")[1] for mergeSamplesListClusters in mergeSamplesListsClusters for cluster in mergeSamplesListClusters]
                
                tsvs = { key : set() for key in mergeSamplesClusters }
                for samplesClusters in self.mergeSamples.values():
                    for cluster in samplesClusters.clusters:
                        clusterNum = cluster.clusterName.split("merged.clusters_")[1]
                        tsvs[clusterNum].add(cluster.tsv)

                mergeClusters = dict.fromkeys(tsvs.keys())
                for clusterNum, tsv in tsvs.items():
                    mergeClusters[clusterNum] = Cluster(clusterName = ("mergedCluster_" + str(clusterNum)), tsv = list(tsvs[clusterNum]), logfile = self.logfile)
                
                self.mergedClusters_results = []
                # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for clusterNum in mergeClusters.keys():
                    outdirConcatenateLanes = os.path.join(outdir, "concatenateLanes")
                    outfile = os.path.join(outdirMergedClusters, ("mergedCluster_" + str(clusterNum) + "_R2.fastq.gz"))
                    
                    clusterFastqs = []
                    for dirpath, subdirs, files in os.walk(outdirConcatenateLanes):
                        for file in files:
                            if file.split("_merged.clusters_")[1].split("_R2.fastq.gz")[0] == clusterNum:
                                clusterFastqs.append(os.path.join(dirpath, file))

                    cmd = clusterFastqs
                    cmd.insert(0, "cat")
                    log.write("Running " + ' '.join(cmd) + "\n\n\n")
                    
                    with open(outfile, "w") as fout:
                        # subprocess.call("date", shell = False, stdout=log, universal_newlines=True)
                        # self.mergedClusters_results.append(executor.submit(subprocess.call, cmd, shell = False, stdout=fout))
                        self.mergedClusters_results.append(subprocess.call(cmd, shell = False, stdout=fout, universal_newlines=True))
                         
                self.mergeSamples = mergeClusters
                self.outdirMergedClusters = outdirMergedClusters
                log.write("self.mergeSamples is now " + str(mergeClusters) + "\n\n\n")
                
                mergedClusters_allSuccess = all(exitCode == 0 for exitCode in self.mergedClusters_results)

                if mergedClusters_allSuccess:
                    msg = "\nmergedClusters finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nmergedClusters did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                    msg = Bcolors.HEADER + "User interrupted. Finishing merging clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                    print(msg)
                    log.write(msg)

    def setMergeClusters(self, outdirMergedClusters):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists(outdirMergedClusters):
                    log.write("Directory not found: " + outdirMergedClusters + " does not exist.")
                    return 2

                mergeSamplesListsClusters = [samplesClusters.clusters for samplesClusters in self.mergeSamples.values()]
                mergeSamplesClusters = [cluster.clusterName.split("merged.clusters_")[1] for mergeSamplesListClusters in mergeSamplesListsClusters for cluster in mergeSamplesListClusters]
                
                tsvs = { key : set() for key in mergeSamplesClusters }
                for samplesClusters in self.mergeSamples.values():
                    for cluster in samplesClusters.clusters:
                        clusterNum = cluster.clusterName.split("merged.clusters_")[1]
                        tsvs[clusterNum].add(cluster.tsv)

                mergeClusters = dict.fromkeys(tsvs.keys())
                for clusterNum, tsv in tsvs.items():
                    mergeClusters[clusterNum] = Cluster(clusterName = ("mergedCluster_" + str(clusterNum)), tsv = list(tsvs[clusterNum]), logfile = self.logfile)
                
                self.mergeSamples = mergeClusters
                self.outdirMergedClusters = outdirMergedClusters
                log.write("self.mergeSamples is now " + str(mergeClusters) + "\n\n\n")
                
            except KeyboardInterrupt:
                    msg = Bcolors.HEADER + "User interrupted. Finishing merging clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                    print(msg)
                    log.write(msg)

    def mapClusters(self, mode, outdir, geneGTF, starIndex, RAM, outTmpDir=None, unique=False, jobs=1):
        print("Running mapClusters with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                self.mapCluster_results = []
                msg = "Mapping clusters.\n"
                log.write(msg)

                if unique:
                    subdirectory = "unique"
                else:
                    subdirectory = "multiple"
                if mode == "merged":
                    samplesDict = self.mergeSamples

                    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                        for clusterNum, cluster in samplesDict.items():
                            fastqdir = os.path.join(outdir, "mergedClusters/")
                            mapOutdir = os.path.join(outdir, "mapCluster/", subdirectory)
                            self.mapCluster_results.append(executor.submit(cluster.mapCluster, "Merged", fastqdir, mapOutdir, geneGTF, starIndex, RAM, outTmpDir, unique, self.slurm, self.modules))
                else:
                    if mode == "perSample":
                        samplesDict = self.samples

                        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                            for sampleId, sample in samplesDict.items():
                                for cluster in sample.clusters:
                                    fastqdir = os.path.join(outdir, "concatenateLanes/", sampleId)
                                    outdir_sample = os.path.join(outdir, "mapCluster/", subdirectory, sampleId)
                                    self.mapCluster_results.append(executor.submit(cluster.mapCluster, sampleId, fastqdir, outdir_sample, geneGTF, starIndex, RAM, outTmpDir, unique, self.slurm, self.modules))
                    else:
                        msg = "Please specify a mode (merged/perSample).\n"
                        print(msg)
                        log.write(msg)
                        return 2

                mapCluster_exitCodes = [i.result()[1] for i in self.mapCluster_results]
                mapCluster_allSuccess = all(exitCode == 0 for exitCode in mapCluster_exitCodes)

                if mapCluster_allSuccess:
                    msg = "\nMapCluster finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nMapCluster did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing mapping for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)
        
    def TEcountsClusters(self, mode, outdir, geneGTF, teGTF, unique=False, jobs=1):
        print("Running TEcounts with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                self.TEcounts_results = []
                msg = "Quantifying TEs.\n"
                log.write(msg)
                
                if unique:
                    subdirectory = "unique"
                else:
                    subdirectory = "multiple"

                if mode == "merged":
                    samplesDict = self.mergeSamples

                    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                        for clusterNum, cluster in samplesDict.items():
                            bam = os.path.join(outdir, "mapCluster/", subdirectory, (cluster.clusterName + "_Aligned.sortedByCoord.out.bam"))
                            outdir_sample = os.path.join(outdir, "TEcounts/", subdirectory)
                            self.TEcounts_results.append(executor.submit(cluster.TEcount, self.name, "Merged", bam, outdir_sample, geneGTF, teGTF, unique, self.slurm, self.modules))
                else:
                    if mode == "perSample":
                        samplesDict = self.samples

                        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                            for sampleId, sample in samplesDict.items():
                                for cluster in sample.clusters:
                                    bam = os.path.join(outdir, "mapCluster/", subdirectory, sampleId, (cluster.clusterName + "_Aligned.sortedByCoord.out.bam"))
                                    outdir_sample = os.path.join(outdir, "TEcounts/", subdirectory, sampleId)
                                    self.TEcounts_results.append(executor.submit(cluster.TEcount, self.name, sampleId, bam, outdir_sample, geneGTF, teGTF, unique, self.slurm, self.modules))
                            
                    else:
                        msg = "Please specify a mode (merged/perSample).\n"
                        print(msg)
                        log.write(msg)
                        return 2
               
                TEcounts_exitCodes = [i.result()[1] for i in self.TEcounts_results]
                TEcounts_allSuccess = all(exitCode == 0 for exitCode in TEcounts_exitCodes)

                if TEcounts_allSuccess:
                    msg = "\nTEcounts finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nTEcounts did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing TEcounts for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def normalizeTECounts(self, mode, outdir, unique=False, jobs=1):
        print("Running normalizeTECounts with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            msg = "Normalizing TE counts.\n"
            log.write(msg)
            try:
                if unique:
                    msg = "\nSorry, normalization only available for multiple mapping.\n"
                    log.write(msg)
                    return False
                else:
                    subdirectory = "multiple"

                indir = os.path.join(outdir, "TEcounts", subdirectory)
                outdirNorm = os.path.join(outdir, "TEcountsNormalized", subdirectory)
                # self.NormalizedOutdir
                
                if mode == "merged":
                    rdata = os.path.join(self.mergeSamplesOutdir, (self.name + ".rds"))

                    if not os.path.exists("mergedSamplesNorm_scripts"):
                        os.makedirs("mergeSamplesNorm_scripts", exist_ok=True)
                    if not os.path.exists(outdirNorm):
                        os.makedirs(outdirNorm, exist_ok=True)
                    
                    cwd = os.path.dirname(os.path.realpath(__file__))
                    cmd = ["Rscript", os.path.join(cwd, "r_scripts/normalize_TEexpression.R"), "-m", mode, "-o", outdirNorm, "-i", indir, "-r", rdata, "-n", self.name]

                    if self.slurmPath != None:
                        cmd = ' '.join(cmd)
                        
                        jobFile =  "mergeSamplesNorm_scripts/" + self.name + "_mergeSamplesNorm.sh"
                        jobId = runJob("normalizeTEcounts", jobFile, cmd, self.slurm, self.modules)
                        
                        msg = sucessSubmit("normalizeTEcounts", self.name, jobId)
                        print(msg)
                        log.write(msg)
                        
                        exitCode = waitForJob(jobId)

                        msg = checkExitCodes("normalizeTEcounts", ("Experiment " + self.name), jobId, exitCode)
                        print(msg)
                        log.write(msg)

                        if exitCode == 0:
                            exitCode = True
                    else:
                        subprocess.call(cmd)
                    self.mergeNormalizedOutdir = outdirNorm

                    if exitCode:
                        msg = "\nTE normalization finished succesfully!\n"
                        log.write(msg)
                        return exitCode
                    else:
                        msg = "\nTE normalization did not finished succesfully.\n"
                        log.write(msg)
                        return exitCode
                else:
                    self.normalized_results = []
                    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                        for sample in list(self.samples.values()):
                            sampleIndir = os.path.join(indir, sample.sampleId)
                            sampleOutdir = os.path.join(outdirNorm, sample.sampleId)
                            self.normalized_results.append(executor.submit(sample.normalizeTEcounts, sampleIndir, sampleOutdir))
                            sample.NormalizedOutdir = sampleOutdir
                    normalized_exitCodes = [i.result() for i in self.normalized_results]
                    normalized_allSuccess = all(exitCode == 0 for exitCode in normalized_exitCodes)

                    if normalized_allSuccess:
                        msg = "\nTE normalization finished succesfully for all samples!\n"
                        log.write(msg)
                        return True
                    else:
                        msg = "\nTE normalization did not finished succesfully for all samples.\n"
                        log.write(msg)
                        return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                print(msg)
                log.write(msg)

    # plotTEexpression -r ../3_mergedSamples/gliomas.RData -m merged -n Gliomas 
    # -t L1HS:L1:LINE,L1PA2:L1:LINE,L1PA3:L1:LINE,L1PA4:L1:LINE,L1PA5:L1:LINE,L1PA6:L1:LINE,L1PA7:L1:LINE,L1PA8:L1:LINE 
    # -i /projects/fs5/raquelgg/Gliomas/Seq073_Seq091/3_mergedSamples/clusterPipeline/TEcountsNormalized 
    # -o /projects/fs5/raquelgg/Gliomas/Seq073_Seq091/3_mergedSamples/clusterPipeline/TEplots
    def plotTEexpression(self, mode, TEsubfamilies, outdir):
        with open(self.logfile, "a") as log:
            msg = "Normalizing TE counts.\n"
            log.write(msg)
            try:
                indir = os.path.join(outdir, "TEcountsNormalized")
                outdirPlots = os.path.join(outdir, "TEplots")
                cwd = os.path.dirname(os.path.realpath(__file__))

                if not os.path.exists("plotTEexpression_scripts"):
                    os.makedirs("plotTEexpression_scripts", exist_ok=True)
                if not os.path.exists(outdirPlots):
                    os.makedirs(outdirPlots, exist_ok=True)

                if mode == "merged":
                    RDatas = os.path.join(self.mergeSamplesOutdir, (self.name + ".rds"))
                    mergedInput = os.path.join(indir, "TE_normalizedValues_aggregatedByClusters_melted.csv")
                    
                    cmd = ["Rscript", os.path.join(cwd, "r_scripts/plot_TEexpression.R"), "-r", RDatas, "-t", TEsubfamilies, "-m", mode, "-n", self.name, "-i", mergedInput, "-o", outdirPlots]

                    if self.slurmPath != None:
                        cmd = ' '.join(cmd)
                        
                        jobFile =  "plotTEexpression_scripts/" + self.name + "_plotTEexpression.sh"
                        jobId = runJob("plotTEexpression", jobFile, cmd, self.slurm, self.modules)
                        
                        msg = sucessSubmit("plotTEexpression", self.name, jobId)
                        print(msg)
                        log.write(msg)
                        
                        exitCode = waitForJob(jobId)

                        msg = checkExitCodes("plotTEexpression", ("Experiment " + self.name), jobId, exitCode)
                        print(msg)
                        log.write(msg)

                        if exitCode == 0:
                            exitCode = True
                    else:
                        subprocess.call(cmd)

                    if exitCode:
                        msg = "\nTE plot finished succesfully!\n"
                        log.write(msg)
                        return exitCode
                    else:
                        msg = "\nTE plot did not finished succesfully.\n"
                        log.write(msg)
                        return exitCode
                else:
                    sampleInputs = ','.join([os.path.join(indir, sample.sampleId, "TE_normalizedValues_melted.csv") for sample in list(self.samples.values())])
                    sampleRDatas = ','.join([sample.rdataPath for sample in list(self.samples.values())])
                    sampleNames = ','.join([sample.sampleId for sample in list(self.samples.values())])
                    modes = ','.join(["perSample" for i in list(self.samples.values())])

                    cmd = ["Rscript", os.path.join(cwd, "r_scripts/plot_TEexpression.R"), "-r", sampleRDatas, "-t", TEsubfamilies, "-m", modes, "-n", sampleNames, "-i", sampleInputs, "-o", outdirPlots, "-c", colourBy]

                    if self.slurmPath != None:
                        cmd = ' '.join(cmd)
                        
                        jobFile =  "plotTEexpression_scripts/perSamples_" + self.name + "_plotTEexpression.sh"
                        jobId = runJob("plotTEexpression", jobFile, cmd, self.slurm, self.modules)
                        
                        msg = sucessSubmit("plotTEexpression", self.name, jobId)
                        print(msg)
                        log.write(msg)
                        
                        exitCode = waitForJob(jobId)

                        msg = checkExitCodes("plotTEexpression", ("Experiment " + self.name), jobId, exitCode)
                        print(msg)
                        log.write(msg)

                        if exitCode == 0:
                            exitCode = True
                    else:
                        subprocess.call(cmd)

                    if exitCode:
                        msg = "\nTE plot finished succesfully!\n"
                        log.write(msg)
                        return exitCode
                    else:
                        msg = "\nTE plot did not finished succesfully.\n"
                        log.write(msg)
                        return exitCode
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                print(msg)
                log.write(msg)

    def processClusters(self, mode, outdir, geneGTF, teGTF, starIndex, RAM, outTmpDir = None, unique=False, jobs=1, finishedTsvToBam = False, finishedFilterUMIs = False, finishedBamToFastq = False, finishedConcatenateLanes = False, finishedMergeClusters = False, finishedMapCluster = False, finishedTEcounts = False, finishedNormalizeTEcounts = False):
        with open(self.logfile, "a") as log:
            msg = "Running whole pipeline.\n"
            log.write(msg)
            finishedOnTheRunMergeClusters = False

            try:
                if mode == "merged":
                    samplesDict = self.mergeSamples
                else:
                    if mode == "perSample":    
                        samplesDict = self.samples
                    else:
                        msg = "Please specify a mode (merged/perSample)"
                        print(msg)
                        log.write(msg)
                        return 2

                if not finishedTsvToBam:
                    current_instruction = "tsvToBam"
                    msg = "Running " + current_instruction
                    log.write(msg)
                    finishedTsvToBam = self.tsvToBamClusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not finishedTsvToBam:
                        msg = "Error in tsvToBam"
                        print(msg)
                        log.write(msg)
                        return 1

                if not finishedFilterUMIs:
                    current_instruction = "filterUMIs"
                    msg = "tsvToBam finished! Moving on to " + current_instruction
                    log.write(msg)
                    finishedFilterUMIs = self.filterUMIsClusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not finishedFilterUMIs:
                        msg = "Error in filterUMIs"
                        print(msg)
                        log.write(msg)
                        return 1
                    
                if not finishedBamToFastq:
                    current_instruction = "bamToFastq"
                    msg = "filterUMIs finished! Moving on to " + current_instruction
                    log.write(msg)
                    finishedBamToFastq = self.bamToFastqClusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not finishedBamToFastq:
                        msg = "Error in BamToFastq"
                        print(msg)
                        log.write(msg)
                        return 1

                if not finishedConcatenateLanes:
                    current_instruction = "concatenateLanes"
                    msg = "bamToFastq finished! Moving on to " + current_instruction
                    log.write(msg)
                    finishedConcatenateLanes = self.concatenateLanesClusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not finishedConcatenateLanes:
                        msg = "Error in ConcatenateLanes"
                        print(msg)
                        log.write(msg)
                        return 1

                if mode == "merged" and not finishedMergeClusters:
                    current_instruction = "mergeClusters"
                    msg = "concatenateLanes finished! You selected merged as your mode, so moving on to " + current_instruction
                    log.write(msg)
                    finishedMergeClusters = self.mergeClusters(outdir = outdir)
                    if not finishedMergeClusters:
                        msg = "Error in mergeClusters"
                        print(msg)
                        log.write(msg)
                        return 1
                    finishedOnTheRunMergeClusters = True

                if not finishedMapCluster:
                    if mode == "merged" and not finishedOnTheRunMergeClusters:
                        self.setMergeClusters(os.path.join(outdir, "mergedClusters"))
                        finishedOnTheRunMergeClusters = True

                    current_instruction = "mapCluster"
                    msg = "mergeClusters finished! Moving on to " + current_instruction
                    log.write(msg)
                    finishedMapCluster = self.mapClusters(mode = mode, outdir = outdir, geneGTF = geneGTF, starIndex = starIndex, RAM = RAM, outTmpDir = outTmpDir, unique = unique, jobs = jobs)
                    if not finishedMapCluster:
                        msg = "Error in MapCluster"
                        print(msg)
                        log.write(msg)
                        return 1

                if not finishedTEcounts:
                    if mode == "merged" and not finishedOnTheRunMergeClusters:
                        self.setMergeClusters(os.path.join(outdir, "mergedClusters"))
                        finishedOnTheRunMergeClusters = True

                    current_instruction = "TEcounts"
                    msg = "mapCluster finished! Moving on to " + current_instruction
                    log.write(msg)
                    finishedTEcounts = self.TEcountsClusters(mode = mode, outdir = outdir, geneGTF = geneGTF, teGTF = teGTF, unique = unique, jobs = jobs)
                    if not finishedTEcounts:
                        msg = "Error in TEcounts"
                        print(msg)
                        log.write(msg)
                        return 1

                if not finishedNormalizeTEcounts:
                    if mode == "merged" and not finishedOnTheRunMergeClusters:
                        self.setMergeClusters(os.path.join(outdir, "mergedClusters"))
                        finishedOnTheRunMergeClusters = True

                    current_instruction = "normalizeTEcounts"
                    msg = "TEcounts finished! Moving on to " + current_instruction
                    log.write(msg)
                    finishedNormalizeTEcounts = self.normalizeTECounts(mode = mode, outdir = outdir, unique = unique, jobs = jobs)
                    if not finishedNormalizeTEcounts:
                        msg = "Error in normalizeTEcounts"
                        print(msg)
                        log.write(msg)
                        return 1

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing instruction " + current_instruction + " for all clusters of all samples before closing." + Bcolors.ENDC + "\n"
                print(msg)
                log.write(msg)

            
