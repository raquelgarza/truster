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
            self.quantifyOutdir = outdir
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + ".\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    def setQuantificationOutdir(self, cellranger_outdir):
        self.quantifyOutdir = cellranger_outdir
        with open(self.logfile, "a") as log:
            msg = "Quantification directory set to: " + cellranger_outdir + ".\n"
            log.write(msg)

    def getClustersAllSamples(self, outdir, excludeFilesPath=None, jobs=1):
        with open(self.logfile, "a") as log:
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample in list(self.samples.values()):
                        sampleIndir = os.path.join(self.quantifyOutdir, sample.sampleId)
                        sampleOutdir = os.path.join(outdir, sample.sampleId)
                        if excludeFilesPath != None:
                            excludeFile = os.path.join(excludeFilesPath, (sample.sampleId + "_exclude.tsv"))
                            if os.path.isfile(excludeFile):
                                msg = "Clustering " + sample.sampleId + " excluding from " + excludeFile + "\n"
                                # args = [sampleIndir, sampleOutdir, excludeFile]
                                # Iterator for more than one argument
                                # executor.submit(lambda p: sample.getClusters(*p), args)
                                executor.submit(sample.getClusters, sampleIndir, sampleOutdir, excludeFile)
                            else:
                                msg = "Clustering " + sample.sampleId + " using all cells. File " + excludeFile + " not found.\n"
                                # args = [sampleIndir, sampleOutdir, None]
                                # executor.submit(lambda p: sample.getClusters(*p), args)
                                executor.submit(sample.getClusters, sampleIndir, sampleOutdir, None)
                            self.getClustersExclusiveOutdir = outdir
                        else:
                            msg = "Clustering " + sample.sampleId + " using all cells.\n"
                            # args = [sampleIndir, sampleOutdir, None]
                            # executor.submit(lambda p: sample.getClusters(*p), args)
                            executor.submit(sample.getClusters, sampleIndir, sampleOutdir, None)
                            self.getClustersOutdir = outdir
                            self.getClustersExclusiveOutdir = None
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
            msg = "The directory for clusters of individual samples is set to: " + clustersOutdir + ".\n"
            log.write(msg)

    def setGetClustersExclusiveOutdir(self, getClustersExclusiveOutdir):
            with open(self.logfile, "a") as log:
                self.getClustersExclusiveOutdir = getClustersExclusiveOutdir
                # Register clusters in each sample 
                for sample in list(self.samples.values()):
                    for i in os.listdir(getClustersExclusiveOutdir):
                        if os.path.isdir(os.path.join(getClustersExclusiveOutdir, i)) and  i == sample.sampleId:
                            msg = sample.registerClustersFromPath(os.path.join(getClustersExclusiveOutdir, sample.sampleId))
                            log.write(msg)
                msg = "The directory for clusters of individual samples is set to: " + getClustersExclusiveOutdir + ". Which includes selected cells. Samples' clustering is now without the excluded cells.\n"
                log.write(msg)

    def velocityAllSamples(self, TE_gtf, geneGTF, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    sampleIndir = os.path.join(self.quantifyOutdir, sample.sampleId)
                    # args = [TE_gtf, geneGTF, sampleIndir]
                    # executor.submit(lambda p: sample.velocity(*p), args)
                    executor.submit(sample.velocity, TE_gtf, geneGTF, sampleIndir)
                executor.shutdown(wait=True)
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + "\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    #def setProcessClustersOutdir(self, processClustersOutdir):
    #    self.processedClustersOutdir = processClustersOutdir

    def mergeSamples(self, outdir):
        # Rscript {input.script} -i {rdata} -n {samplenames} -o {params.outpath}
        # Paths to RData files
        with open(self.logfile, "a") as log:
            if self.getClustersExclusiveOutdir != None:
                samplesSeuratRdata = [os.path.join(self.getClustersExclusiveOutdir, sample.sampleId, (sample.sampleId + ".RData")) for sample in list(self.samples.values())]
            else:
                samplesSeuratRdata = [os.path.join(self.getClustersOutdir, sample.sampleId, (sample.sampleId + ".RData")) for sample in list(self.samples.values())]
            
            # Sample ids
            samplesIds = [sample.sampleId for sample in list(self.samples.values())]
            
            msg = "Merging samples to produce a combined clustering.\n"
            log.write(msg)
            # If we haven't made the merge before, create a directory to store the scripts needed to do so
            if not os.path.exists("mergeSamples_scripts"):
                os.makedirs("mergeSamples_scripts", exist_ok=True)
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
                
            # Run script ../r_scripts/merge_samples.R with input (-i) of the RData paths
            # and output (-o) of the output directory desired, -s for sample ids,
            # and -e for sample names used in cellranger
            cwd = os.path.dirname(os.path.realpath(__file__))
            cmd = ["Rscript", os.path.join(cwd, "r_scripts/merge_samples.R"), "-i", ','.join(samplesSeuratRdata), "-o", outdir, "-s", ','.join(samplesIds), "-e", self.name]
            
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
                            bam = os.path.join(self.quantifyOutdir, sampleId, "outs/possorted_genome_bam.bam")
                            outdir_sample = os.path.join(outdir, "tsvToBam/", sampleId)
                            # args = [sampleId, bam, outdir_sample, self.slurm, self.modules]
                            # self.tsvToBam_results.append(executor.submit(lambda p: cluster.tsvToBam(*p), args))
                            exitCode = executor.submit(cluster.tsvToBam, sampleId, bam, outdir_sample, self.slurm, self.modules)
                            self.tsvToBam_results.append(exitCode)
                            msg = cluster.clusterName + " return " + str(exitCode.result()[1]) + " from tsvToBam.\n"
                            log.write(msg)
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
                            #clusters = self.mergeSamplesClusters[sampleId]
                            # args = [sampleId, inbam, outdir_sample, self.slurm, self.modules]
                            # self.filterUMIs_results.append(executor.submit(lambda p: cluster.filterUMIs(*p), args))
                            exitCode = executor.submit(cluster.filterUMIs, sampleId, inbam, outdir_sample, self.slurm, self.modules)
                            self.filterUMIs_results.append(exitCode)
                            msg = cluster.clusterName + " return " + str(exitCode.result()[1]) + " from filterUMIs.\n"
                            log.write(msg)
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
                            # args = [sampleId, bam, outdir_sample, self.slurm, self.modules]
                            # self.bamToFastq_results.append(executor.submit(lambda p: cluster.bamToFastq(*p), args))
                            exitCode = executor.submit(cluster.bamToFastq, sampleId, bam, outdir_sample, self.slurm, self.modules)
                            self.bamToFastq_results.append(exitCode)
                            msg = cluster.clusterName + " return " + str(exitCode.result()[1]) + " from bamToFastq.\n"
                            log.write(msg)
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
                            # args = [sampleId, indir, outdir_sample, self.slurm, self.modules]
                            # self.concatenateLanes_results.append(executor.submit(lambda p: cluster.concatenateLanes(*p), args))
                            exitCode = executor.submit(cluster.concatenateLanes, sampleId, indir, outdir_sample, self.slurm, self.modules)
                            self.concatenateLanes_results.append(exitCode)
                            msg = cluster.clusterName + " return " + str(exitCode.result()[1]) + " from concatenateLanes.\n"
                            log.write(msg)
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

    def mapClusters(self, mode, outdir, geneGTF, starIndex, jobs=1):
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

                self.mapCluster_results = []
                msg = "Mapping clusters.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            fastqdir = os.path.join(outdir, "concatenateLanes/", sampleId)
                            outdir_sample = os.path.join(outdir, "mapCluster/", sampleId)
                            # args = [sampleId, fastqdir, outdir_sample, geneGTF, starIndex, self.slurm, self.modules]
                            # self.mapCluster_results.append(executor.submit(lambda p: cluster.mapCluster(*p), args))
                            exitCode = executor.submit(cluster.mapCluster, sampleId, fastqdir, outdir_sample, geneGTF, starIndex, self.slurm, self.modules)
                            self.mapCluster_results.append(exitCode)
                            msg = cluster.clusterName + " return " + str(exitCode.result()[1]) + " from mapCluster.\n"
                            log.write(msg)
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
        
    def TEcountsClusters(self, mode, outdir, geneGTF, teGTF, jobs=1):
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

                self.TEcounts_results = []
                msg = "Quantifying TEs.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            bam = os.path.join(outdir, "mapCluster/", sampleId, (cluster.clusterName + "_Aligned.sortedByCoord.out.bam"))
                            outdir_sample = os.path.join(outdir, "TEcounts/", sampleId)
                            # args = [self.name, sampleId, bam, outdir_sample, geneGTF, teGTF, self.slurm, self.modules]
                            # self.TEcounts_results.append(executor.submit(lambda p: cluster.TEcount(*p), args))
                            exitCode = executor.submit(cluster.TEcount, self.name, sampleId, bam, outdir_sample, geneGTF, teGTF, self.slurm, self.modules)
                            self.TEcounts_results.append(exitCode)
                            msg = cluster.clusterName + " return " + str(exitCode.result()[1]) + " from TEcount.\n"
                            log.write(msg)
                TEcounts_exitCodes = [i.result()[1] for i in self.TEcounts_results]
                TEcounts_allSuccess = all(exitCode == 0 for exitCode in TEcounts_exitCodes)

                if TEcountsCluster_allSuccess:
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

    def normalizeTECountsCluster(self, mode, jobs=1):
        with open(self.logfile, "a") as log:
            msg = "Normalizing TE counts.\n"
            log.write(msg)
            try:
                if mode == "merged":
                    rdata = os.path.join(self.mergeSamplesOutdir, (self.name + ".RData"))
                    indir = os.path.join(self.mergeSamplesOutdir, "clusterPipeline/TEcounts")
                    outdir = os.path.join(self.mergeSamplesOutdir, "clusterPipeline/TEcountsNormalized")
                    self.NormalizedOutdir = outdir
                     
                    if not os.path.exists("mergedSamplesNorm_scripts"):
                        os.makedirs("mergeSamplesNorm_scripts", exist_ok=True)
                    if not os.path.exists(outdir):
                        os.makedirs(outdir, exist_ok=True)
                    
                    cwd = os.path.dirname(os.path.realpath(__file__))
                    cmd = ["Rscript", os.path.join("r_scripts/normalize_TEexpression.R"), "-m", mode, "-o", outdir, "-i", indir, "-r", rdata]

                    if self.slurmPath != None:
                        cmd = ' '.join(cmd)
                        
                        jobFile =  "mergeSamplesNorm_scripts/" + self.name + "_mergeSamplesNorm.sh"
                        
                        jobId = runJob("mergeSamplesNorm", jobFile, cmd, self.slurm, self.modules)
                        
                        msg = sucessSubmit("mergeSamplesNorm", self.name, jobId)
                        print(msg)
                        log.write(msg)
                        
                        exitCode = waitForJob(jobId)
                        msg = checkExitCodes("mergeSamplesNorm", ("Experiment " + self.name), jobId, exitCode)
                        print(msg)
                        log.write(msg)

                    else:
                        subprocess.call(cmd)
                    self.mergeNormalizedTEcountsOutdir = outdir

                    if exitCode:
                        msg = "\nTE normalization finished succesfully!\n"
                        log.write(msg)
                        return True
                    else:
                        msg = "\nTE normalization did not finished succesfully.\n"
                        log.write(msg)
                        return False
                else:
                    self.normalized_results = []
                    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                        for sample in list(self.samples.values()):
                            sampleIndir = os.path.join(indir, sample.sampleId)
                            sampleOutdir = os.path.join(outdir, sample.sampleId)
                            # args = [sampleIndir, sampleOutdir]
                            # self.normalized_results.append(executor.submit(lambda p: sample.normalizeTEcounts(*p), args))
                            exitCode = executor.submit(sample.normalizeTEcounts, sampleIndir, sampleOutdir)
                            self.normalized_results.append(exitCode)
                            msg = cluster.clusterName + " return " + str(exitCode.result()[1]) + " from TE count normalization.\n"
                            log.write(msg)
                            sample.NormalizedOutdir = outdir
                    normalized_exitCodes = [i.result()[1] for i in self.normalized_results]
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

    def processClusters(self, mode, outdir, geneGTF, teGTF, starIndex, jobs=1):
        with open(self.logfile, "a") as log:
            msg = "Running whole pipeline.\n"
            log.write(msg)
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

                current_instruction = "tsvToBam"
                msg = "Running " + current_instruction
                log.write(msg)
                
                # tsvToBam_allSuccess = self.tsvToBamClusters(mode, outdir, jobs)
                tsvToBam_allSuccess = True

                if tsvToBam_allSuccess:
                    current_instruction = "filterUMIs"
                    msg = "tsvToBam finished! Moving on to " + current_instruction
                    log.write(msg)
                    filterUMIs_allSuccess = self.filterUMIsClusters(mode, outdir, jobs)

                    if filterUMIs_allSuccess:
                        current_instruction = "bamToFastq"
                        msg = "filterUMIs finished! Moving on to " + current_instruction
                        log.write(msg)
                        bamToFastq_allSuccess = self.bamToFastqClusters(mode, outdir, jobs)
        
                        if bamToFastq_allSuccess:
                            current_instruction = "concatenateLanes"
                            msg = "bamToFastq finished! Moving on to " + current_instruction
                            log.write(msg)
                            concatenateLanes_allSuccess = self.concatenateLanesClusters(mode, outdir, jobs)

                            # self.concatenateLanes_results = []
                            # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                            #     for sampleId, sample in samplesDict.items():
                            #         for cluster in sample.clusters:
                            #             indir = os.path.join(outdir, "bamToFastq/", sampleId)
                            #             outdir_sample = os.path.join(outdir, "concatenateLanes/", sampleId)
                            #             self.concatenateLanes_results.append(executor.submit(cluster.concatenateLanes, sampleId, indir, outdir_sample, self.slurm, self.modules))
                            # concatenateLanes_exitCodes = [i.result()[1] for i in self.concatenateLanes_results]
                            # concatenateLanes_allSuccess = all(exitCode == 0 for exitCode in concatenateLanes_exitCodes)

                            if concatenateLanes_allSuccess:
                                current_instruction = "mapCluster"
                                msg = "concatenateLanes finished! Moving on to " + current_instruction
                                log.write(msg)
                                mapCluster_allSuccess = self.mapClusters(mode, outdir, geneGTF, starIndex, jobs)

                                # self.mapCluster_results = []
                                # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                                #     for sampleId, sample in samplesDict.items():
                                #         for cluster in sample.clusters:
                                #             fastqdir = os.path.join(outdir, "concatenateLanes/", sampleId)
                                #             outdir_sample = os.path.join(outdir, "mapCluster/", sampleId)
                                #             self.mapCluster_results.append(executor.submit(cluster.mapCluster, sampleId, fastqdir, outdir_sample, geneGTF, starIndex, self.slurm, self.modules))
                                # mapCluster_exitCodes = [i.result()[1] for i in self.mapCluster_results]
                                # mapCluster_allSuccess = all(exitCode == 0 for exitCode in mapCluster_exitCodes)
                                if mapCluster_allSuccess:
                                    current_instruction = "TEcounts"
                                    msg = "mapCluster finished! Moving on to " + current_instruction
                                    log.write(msg)
                                    TEcounts_allSuccess = self.TEcountsClusters(mode, outdir, geneGTF, teGTF, jobs)

                                    # self.TEcounts_results = []
                                    # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                                    #     for sampleId, sample in samplesDict.items():
                                    #         for cluster in sample.clusters:
                                    #             bam = os.path.join(outdir, "mapCluster/", sampleId, (cluster.clusterName + "_Aligned.sortedByCoord.out.bam"))
                                    #             outdir_sample = os.path.join(outdir, "TEcounts/", sampleId)
                                    #             self.TEcounts_results.append(executor.submit(cluster.TEcount, self.name, sampleId, bam, outdir_sample, geneGTF, teGTF, self.slurm, self.modules))
                                    # TEcounts_exitCodes = [i.result()[1] for i in self.TEcounts_results]
                                    # TEcounts_allSuccess = all(exitCode == 0 for exitCode in TEcounts_exitCodes)
                                    if TEcounts_allSuccess:
                                        current_instruction = "normalizeTEcounts"
                                        msg = "TEcounts finished! Moving on to " + current_instruction
                                        log.write(msg)
                                        normalizeTEcounts_allSuccess = self.normalizeTECountsCluster(mode, jobs)
                                    else:
                                        msg = "Error in normalizeTEcounts"
                                        print(msg)
                                        log.write(msg)
                                        return 1
                                else:
                                    msg = "Error in mapCluster"
                                    print(msg)
                                    log.write(msg)
                                    return 1
                            else:
                                msg = "Error in concatenateLanes"
                                print(msg)
                                log.write(msg)
                                return 1
                        else:
                            msg = "Error in bamToFastq"
                            print(msg)
                            log.write(msg)
                            return 1
                    else:
                        msg = "Error in filterUMIs"
                        print(msg)
                        log.write(msg)
                        return 1
                else:
                    msg = "Error in tsvToBam"
                    print(msg)
                    log.write(msg)
                    return 1
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing instruction " + current_instruction + " for all clusters of all samples before closing." + Bcolors.ENDC + "\n"
                print(msg)
                log.write(msg)

        
