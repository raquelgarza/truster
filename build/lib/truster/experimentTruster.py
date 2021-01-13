import subprocess
import json
from .jobHandler import *
import os
from .sampleTruster import sampleTruster
from .clusterTruster import clusterTruster
from .bcolors import bcolors
import concurrent.futures
import copy

class experimentTruster:

    def __init__(self, name="", slurmPath=None, modulesPath=None):
        self.name = name
        self.slurmPath = slurmPath
        self.modulesPath = modulesPath
        self.samples = {}

        if slurmPath != None:
            try:
                with open(slurmPath, "r") as config_file:
                    self.slurm = json.load(config_file)
            except FileNotFoundError:
                print(bcolors.FAIL + "Error: Cluster configuration file not found")
                sys.exit(1)

        if modulesPath != None:
            try:
                with open(modulesPath, "r") as modules_file:
                    self.modules = json.load(modules_file)
            except FileNotFoundError:
                print(bcolors.FAIL + "Error: Module configuration file not found")
                sys.exit(1)

    def registerSample(self, sampleId = "", sampleName = "", rawPath = ""):
        newSample = {sampleId:sampleTruster(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = rawPath)}
        self.samples = {**self.samples, **newSample}
        # self.samples.extend(sampleTruster(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName))

    def unregisterSample(self, sampleId):
        self.samples.pop(sampleId)

    def registerSamplesFromPath(self, indir=[], folderNamesAsSampleIds=True):
        if(folderNamesAsSampleIds):
            for i in indir:
                samples = {item[0] : item[1].split("_L00")[0] for item in [(os.path.join(dp, f)).split("/")[-2:] for dp, dn, fn in os.walk(os.path.expanduser(i)) for f in fn]}
                newSamples = {sampleId:sampleTruster(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = i) for sampleId, sampleName in samples.items()}
                self.samples = {**self.samples, **newSamples}
                # self.samples.extend([sampleTruster(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = i) for sampleId, sampleName in samples.items()])

    def quantify(self, crIndex, outdir, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:

                for sample in self.samples.values():
                    sampleIndir = os.path.join(sample.rawPath, sample.sampleId)
                    sampleOutdir = os.path.join(outdir, sample.sampleId)

                    executor.submit(sample.quantify, crIndex, sampleIndir, sampleOutdir)
            self.quantifyOutdir = outdir
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

    def setQuantificationOutdir(self, cellranger_outdir):
        self.quantifyOutdir = cellranger_outdir

    def getClustersAllSamples(self, outdir, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    sampleIndir = os.path.join(self.quantifyOutdir, sample.sampleId)
                    sampleOutdir = os.path.join(outdir, sample.sampleId)
                    executor.submit(sample.getClusters, sampleIndir, sampleOutdir)
                executor.shutdown(wait=True)
            self.getClustersOutdir = outdir
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

    def setClusteringOutdir(self, seuratRdata_outdir):
        self.getClustersOutdir = seuratRdata_outdir
        # Add cluster objects

    def velocityAllSamples(self, TE_gtf, geneGTF, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    sampleIndir = os.path.join(self.quantifyOutdir, sample.sampleId)
                    executor.submit(sample.velocity, TE_gtf, geneGTF, sampleIndir)
                executor.shutdown(wait=True)
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

    def setProcessClustersOutdir(self, processClustersOutdir):
        self.processedClustersOutdir = processClustersOutdir

    def mergeSamples(self, outdir):
        # Rscript {input.script} -i {rdata} -n {samplenames} -o {params.outpath}
        # Paths to RData files
        samplesSeuratRdata = [os.path.join(self.getClustersOutdir, sample.sampleId, (sample.sampleId + ".RData")) for sample in list(self.samples.values())]
        
            # Sample ids
        samplesIds = [sample.sampleId for sample in list(self.samples.values())]
        
            # If we haven't made the merge before, create a directory to store the scripts needed to do so
        if not os.path.exists("mergeSamples_scripts"):
            os.makedirs("mergeSamples_scripts", exist_ok=True)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
            
        # Run script ../r_scripts/merge_samples.R with input (-i) of the RData paths
        # and output (-o) of the output directory desired, -s for sample ids,
        # and -e for sample names used in cellranger
        cmd = ["Rscript", "../r_scripts/merge_samples.R", "-i", ','.join(samplesSeuratRdata), "-o", outdir, "-s", ','.join(samplesIds), "-e", self.name]
        
        # If we are on a server with slurm, use the configuration file to send the jobs 
        if self.slurmPath != None:    
            cmd = ' '.join(cmd)
            
            # This command will be stored along with the slurm configurations at
            # mergeSamples_scripts/sample_mergeSamples.sh
            jobFile =  "mergeSamples_scripts/" + self.name + "_mergeSamples.sh"
            
            try:
                # Run job script
                jobId = jobHandler.runJob("mergeSamples", jobFile, cmd, self.slurm, self.modules)
                
                           # Print if the job was succesfully submitted
                print(jobHandler.sucessSubmit("mergeSamples", self.name, jobId))
                
                        # Wait for the job to finish and returns an exit code
                exitCode = jobHandler.waitForJob(jobId)
            
                # Print if the job was finished succesfully or not
                print(jobHandler.checkExitCodes("mergeSamples", ("Experiment " + self.name), jobId, exitCode))
                
                       # If it finished succesfully then 
                if exitCode == 0:
            
                               # For each of the registered samples
                    self.mergeSamples = copy.deepcopy(self.samples)

                    # Empty the clusters bc we made new ones (shared/merged)
                    for k,v in self.mergeSamples.items():
                        v.emptyClusters()
                
                    print("Emptied clusters")
                
                    print([j.clusters for j in self.mergeSamples.values()])

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
                            cluster = clusterTruster(clusterName, os.path.join(outdir, i))
                            print(sampleId, cluster)
                            print(self.mergeSamples[sampleId].clusters)
                            self.mergeSamples[sampleId].clusters.append(cluster)
                    print([j.clusters for j in self.mergeSamples.values()])
                return exitCode
            except:
                print(jobHandler.genericError("mergeSamples", self.name))
                return
        else:
            subprocess.call(cmd)
            
        self.mergeSamplesOutdir = outdir

    def tsvToBamClusters(self, mode, outdir, jobs=1):
        try:
            if mode == "merged":
                samplesDict = self.mergeSamples
            else:
                if mode == "perSample":
                    samplesDict = self.samples
                else:
                    print("Please specify a mode (merged/perSample)")
                    return
                
            self.tsvToBam_results = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sampleId, sample in samplesDict.items():
                    for cluster in sample.clusters:
                        bam = os.path.join(self.quantifyOutdir, sampleId, "outs/possorted_genome_bam.bam")
                        outdir_sample = os.path.join(outdir, "tsvToBam/", sampleId)
                        self.tsvToBam_results.append(executor.submit(cluster.tsvToBam, sampleId, bam, outdir_sample, self.slurm, self.modules))
            tsvToBam_exitCodes = [i.result()[1] for i in self.tsvToBam_results]
            tsvToBam_allSuccess = all(exitCode == 0 for exitCode in tsvToBam_exitCodes)

            if tsvToBam_allSuccess:
                return 0
            else:
                return tsvToBam_exitCodes

        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted. Finishing tsvToBam for all clusters of all samples before closing." + bcolors.ENDC)
 

    def filterUMIsClusters(self, mode, outdir, jobs=1):
        try:
            if mode == "merged":
                samplesDict = self.mergeSamples
            else:
                if mode == "perSample":
                    samplesDict = self.samples
                else:
                    print("Please specify a mode (merged/perSample)")
                    return

                self.filterUMIs_results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            inbam = os.path.join(outdir, "tsvToBam/", sampleId, (cluster.clusterName + ".bam"))
                            outdir_sample = os.path.join(outdir, "filterUMIs/", sampleId)
                            #clusters = self.mergeSamplesClusters[sampleId]
                            self.filterUMIs_results.append(executor.submit(cluster.filterUMIs, sampleId, inbam, outdir_sample, self.slurm, self.modules))
                filterUMIs_exitCodes = [i.result()[1] for i in self.filterUMIs_results]
                filterUMIs_allSuccess = all(exitCode == 0 for exitCode in filterUMIs_exitCodes)

            if filterUMIs_allSuccess:
                return 0
            else:
                return filterUMIs_exitCodes
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted. Finishing filterUMIs for all clusters of all samples before closing." + bcolors.ENDC)
 
    def bamToFastqClusters(self, mode, outdir, jobs=1):
        try:
            if mode == "merged":
                samplesDict = self.mergeSamples
            else:
                if mode == "perSample":
                    samplesDict = self.samples
                else:
                    print("Please specify a mode (merged/perSample)")
                    return
            self.bamToFastq_results = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sampleId, sample in samplesDict.items():
                    for cluster in sample.clusters:
                        bam = os.path.join(outdir, "filterUMIs/", sampleId, (cluster.clusterName + "_filtered.bam"))
                        outdir_sample = os.path.join(outdir, "bamToFastq/", sampleId)
                        self.bamToFastq_results.append(executor.submit(cluster.bamToFastq, sampleId, bam, outdir_sample, self.slurm, self.modules))
            bamToFastq_exitCodes = [i.result()[1] for i in self.bamToFastq_results]
            bamToFastq_allSuccess = all(exitCode == 0 for exitCode in bamToFastq_exitCodes)

            if bamToFastq_allSuccess:
                return 0
            else:
                return bamToFastq_exitCodes
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted. Finishing bamToFastq for all clusters of all samples before closing." + bcolors.ENDC)

    def concatenateLanesClusters(self, mode, outdir, jobs=1):
        try:
            if mode == "merged":
                samplesDict = self.mergeSamples
            else:
                if mode == "perSample":
                    samplesDict = self.samples
                else:
                    print("Please specify a mode (merged/perSample)")
                    return

            self.concatenateLanes_results = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sampleId, sample in samplesDict.items():
                    for cluster in sample.clusters:
                        indir = os.path.join(outdir, "bamToFastq/", sampleId)
                        outdir_sample = os.path.join(outdir, "concatenateLanes/", sampleId)
                        self.concatenateLanes_results.append(executor.submit(cluster.concatenateLanes, sampleId, indir, outdir_sample, self.slurm, self.modules))
            concatenateLanes_exitCodes = [i.result()[1] for i in self.concatenateLanes_results]
            concatenateLanes_allSuccess = all(exitCode == 0 for exitCode in concatenateLanes_exitCodes)

            if concatenateLanes_allSuccess:
                return 0
            else:
                return concatenateLanes_exitCodes
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted. Finishing concatenateLanes for all clusters of all samples before closing." + bcolors.ENDC)

    def mapClusters(self, mode, outdir, geneGTF, starIndex, jobs=1):
        try:
            if mode == "merged":
                samplesDict = self.mergeSamples
            else:
                if mode == "perSample":
                    samplesDict = self.samples
                else:
                    print("Please specify a mode (merged/perSample)")
                    return

            self.mapCluster_results = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sampleId, sample in samplesDict.items():
                    for cluster in sample.clusters:
                        fastqdir = os.path.join(outdir, "concatenateLanes/", sampleId)
                        outdir_sample = os.path.join(outdir, "mapCluster/", sampleId)
                        self.mapCluster_results.append(executor.submit(cluster.mapCluster, sampleId, fastqdir, outdir_sample, geneGTF, starIndex, self.slurm, self.modules))
            mapCluster_exitCodes = [i.result()[1] for i in self.mapCluster_results]
            mapCluster_allSuccess = all(exitCode == 0 for exitCode in mapCluster_exitCodes)

            if mapCluster_allSuccess:
                return 0
            else:
                return mapCluster_exitCodes
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted. Finishing mapping for all clusters of all samples before closing." + bcolors.ENDC)
        

    # TEcount(self, experimentName, sampleId, bam, outdir, geneGTF, teGTF, slurm=None, modules=None):
    def TEcountsClusters(self, mode, outdir, geneGTF, teGTF, jobs=1):
        try:
            if mode == "merged":
                samplesDict = self.mergeSamples
            else:
                if mode == "perSample":
                    samplesDict = self.samples
                else:
                    print("Please specify a mode (merged/perSample)")
                    return

            self.TEcounts_results = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sampleId, sample in samplesDict.items():
                    for cluster in sample.clusters:
                        bam = os.path.join(outdir, "mapCluster/", sampleId, (cluster.clusterName + "_Aligned.sortedByCoord.out.bam"))
                        outdir_sample = os.path.join(outdir, "TEcounts/", sampleId)
                        self.TEcounts_results.append(executor.submit(cluster.TEcount, self.name, sampleId, bam, outdir_sample, geneGTF, teGTF, self.slurm, self.modules))
            TEcounts_exitCodes = [i.result()[1] for i in self.TEcounts_results]
            TEcounts_allSuccess = all(exitCode == 0 for exitCode in TEcounts_exitCodes)

            if TEcountsCluster_allSuccess:
                return 0
            else:
                return TEcountsCluster_exitCodes

        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted. Finishing TEcounts for all clusters of all samples before closing." + bcolors.ENDC)
        

    def processClusters(self, mode, outdir, geneGTF, teGTF, starIndex, jobs=1):
        try:
            if mode == "merged":
                samplesDict = self.mergeSamples
            else:
                if mode == "perSample":    
                        samplesDict = self.samples
                else:
                    print("Please specify a mode (merged/perSample)")
                    return

            current_instruction = "tsvToBam"
            # tsvToBam(self, sampleId, bam, outdir, slurm=None, modules=None):
            self.tsvToBam_results = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sampleId, sample in samplesDict.items():
                    for cluster in sample.clusters:
                        bam = os.path.join(self.quantifyOutdir, sampleId, "outs/possorted_genome_bam.bam")
                        outdir_sample = os.path.join(outdir, "tsvToBam/", sampleId)
                        self.tsvToBam_results.append(executor.submit(cluster.tsvToBam, sampleId, bam, outdir_sample, self.slurm, self.modules))
            tsvToBam_exitCodes = [i.result()[1] for i in self.tsvToBam_results]
            tsvToBam_allSuccess = all(exitCode == 0 for exitCode in tsvToBam_exitCodes)

            if tsvToBam_allSuccess:
                current_instruction = "filterUMIs"
                # filterUMIs(self, sampleId, inbam, outdir, slurm=None, modules=None)
                self.filterUMIs_results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sampleId, sample in samplesDict.items():
                        for cluster in sample.clusters:
                            inbam = os.path.join(outdir, "tsvToBam/", sampleId, (cluster.clusterName + ".bam"))
                            outdir_sample = os.path.join(outdir, "filterUMIs/", sampleId)
                            #clusters = self.mergeSamplesClusters[sampleId]
                            self.filterUMIs_results.append(executor.submit(cluster.filterUMIs, sampleId, inbam, outdir_sample, self.slurm, self.modules))
                filterUMIs_exitCodes = [i.result()[1] for i in self.filterUMIs_results]
                filterUMIs_allSuccess = all(exitCode == 0 for exitCode in filterUMIs_exitCodes)
                if filterUMIs_allSuccess:
                    current_instruction = "bamToFastq"
                    # bamToFastq(self, sampleId, bam, outdir, slurm=None, modules=None):
                    self.bamToFastq_results = []
                    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                        for sampleId, sample in samplesDict.items():
                            for cluster in sample.clusters:
                                bam = os.path.join(outdir, "filterUMIs/", sampleId, (cluster.clusterName + "_filtered.bam"))
                                outdir_sample = os.path.join(outdir, "bamToFastq/", sampleId)
                                self.bamToFastq_results.append(executor.submit(cluster.bamToFastq, sampleId, bam, outdir_sample, self.slurm, self.modules))
                    bamToFastq_exitCodes = [i.result()[1] for i in self.bamToFastq_results]
                    bamToFastq_allSuccess = all(exitCode == 0 for exitCode in bamToFastq_exitCodes)
    
                    if bamToFastq_allSuccess:
                        current_instruction = "concatenateLanes"
                        # concatenateLanes(self, sampleId, indir, outdir, slurm=None, modules=None):
                        self.concatenateLanes_results = []
                        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                            for sampleId, sample in samplesDict.items():
                                for cluster in sample.clusters:
                                    indir = os.path.join(outdir, "bamToFastq/", sampleId)
                                    outdir_sample = os.path.join(outdir, "concatenateLanes/", sampleId)
                                    self.concatenateLanes_results.append(executor.submit(cluster.concatenateLanes, sampleId, indir, outdir_sample, self.slurm, self.modules))
                        concatenateLanes_exitCodes = [i.result()[1] for i in self.concatenateLanes_results]
                        concatenateLanes_allSuccess = all(exitCode == 0 for exitCode in concatenateLanes_exitCodes)

                        if concatenateLanes_allSuccess:
                            current_instruction = "mapCluster"
                            # mapCluster(self, sampleId, fastqdir, outdir, geneGTF, starIndex, slurm=None, modules=None):
                            self.mapCluster_results = []
                            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                                for sampleId, sample in samplesDict.items():
                                    for cluster in sample.clusters:
                                        fastqdir = os.path.join(outdir, "concatenateLanes/", sampleId)
                                        outdir_sample = os.path.join(outdir, "mapCluster/", sampleId)
                                        self.mapCluster_results.append(executor.submit(cluster.mapCluster, sampleId, fastqdir, outdir_sample, geneGTF, starIndex, self.slurm, self.modules))
                            mapCluster_exitCodes = [i.result()[1] for i in self.mapCluster_results]
                            mapCluster_allSuccess = all(exitCode == 0 for exitCode in mapCluster_exitCodes)
                            if mapCluster_allSuccess:
                                current_instruction = "TEcounts"
                                # TEcount(self, experimentName, sampleId, bam, outdir, geneGTF, teGTF, slurm=None, modules=None):
                                self.TEcounts_results = []
                                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                                    for sampleId, sample in samplesDict.items():
                                        for cluster in sample.clusters:
                                            bam = os.path.join(outdir, "mapCluster/", sampleId, (cluster.clusterName + "_Aligned.sortedByCoord.out.bam"))
                                            outdir_sample = os.path.join(outdir, "TEcounts/", sampleId)
                                            self.TEcounts_results.append(executor.submit(cluster.TEcount, self.name, sampleId, bam, outdir_sample, geneGTF, teGTF, self.slurm, self.modules))
                                TEcounts_exitCodes = [i.result()[1] for i in self.TEcounts_results]
                                TEcounts_allSuccess = all(exitCode == 0 for exitCode in TEcounts_exitCodes)
                            else:
                                print("Error in mapCluster")
                                return 1
                        else:
                            print("Error in concatenateLanes")
                            return 1
                    else:
                        print("Error in bamToFastq")
                        return 1
                else:
                    print("Error in filterUMIs")
                    return 1
            else:
                print("Error in tsvToBam")
                return 1

        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted. Finishing instruction " + current_instruction + " for all clusters of all samples before closing." + bcolors.ENDC)
        

        if mode == "merged":
            self.processedClustersMergedSamplesOutdir = outdir
        else:
            self.processedClustersOutdir = outdir


    def setProcessedClustersMergedSamplesOutdir(self, processedClustersMergedSamplesOutdir):
        self.processedClustersMergedSamplesOutdir = processedClustersMergedSamplesOutdi

    def setMergeSamplesOutdir(self, mergeSamplesOutdir):
        self.mergeSamplesOutdir = mergeSamplesOutdir


    def normalizeTECounts(self, mode, jobs=1):
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
                
                cmd = ["Rscript", "../r_scripts/normalize_TEexpression.R", "-m", mode, "-o", outdir, "-i", indir, "-r", rdata]

                if self.slurmPath != None:
                    cmd = ' '.join(cmd)
                    
                    jobFile =  "mergeSamplesNorm_scripts/" + self.name + "_mergeSamplesNorm.sh"
                    
                    jobId = jobHandler.runJob("mergeSamplesNorm", jobFile, cmd, self.slurm, self.modules)
                    
                    print(jobHandler.sucessSubmit("mergeSamplesNorm", self.name, jobId))
                    
                    exitCode = jobHandler.waitForJob(jobId)
                    print(jobHandler.checkExitCodes("mergeSamplesNorm", ("Experiment " + self.name), jobId, exitCode))
                else:
                    subprocess.call(cmd)
                self.mergeNormalizedTEcountsOutdir = outdir

            else:
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample in list(self.samples.values()):
                        sampleIndir = os.path.join(indir, sample.sampleId)
                        sampleOutdir = os.path.join(outdir, sample.sampleId)
                        executor.submit(sample.normalizeTEcounts, sampleIndir, sampleOutdir)
                        sample.NormalizedOutdir = outdir
                    executor.shutdown(wait=True)

        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)


