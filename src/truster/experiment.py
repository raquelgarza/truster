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

        with open(self.logfile, "a") as log:
	        if slurmPath != None:
	            try:
	                with open(slurmPath, "r") as config_file:
	                    self.slurm = json.load(config_file)
	            except FileNotFoundError:
	                msg = bcolors.FAIL + "Error: Cluster configuration file not found" + bcolors.ENDC
					print(msg)
					log.write(msg)
	                return 1

	        if modulesPath != None:
            try:
                with open(modulesPath, "r") as modules_file:
                    self.modules = json.load(modules_file)
            except FileNotFoundError:
                msg = bcolors.FAIL + "Error: Module configuration file not found" + bcolors.ENDC
				print(msg)
				log.write(msg)
                sys.exit(1)

    def registerSample(self, sampleId = "", sampleName = "", rawPath = ""):
        newSample = {sampleId : Sample(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = rawPath)}
        self.samples = {**self.samples, **newSample}
        # self.samples.extend(sampleTruster(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName))

    def unregisterSample(self, sampleId):
        self.samples.pop(sampleId)

    def registerSamplesFromPath(self, indir=[], folderNamesAsSampleIds=True):
        if(folderNamesAsSampleIds):
            for i in indir:
                samples = {item[0] : item[1].split("_L00")[0] for item in [(os.path.join(dp, f)).split("/")[-2:] for dp, dn, fn in os.walk(os.path.expanduser(i)) for f in fn]}
                newSamples = {sampleId : Sample(slurm=self.slurm, modules = self.modules, sampleId = sampleId, sampleName = sampleName, rawPath = i) for sampleId, sampleName in samples.items()}
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
            msg = bcolors.HEADER + "User interrupted" + bcolors.ENDC
            print(msg)
            with open(self.logfile, "a") as log:
            	log.write(msg)

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
            msg = bcolors.HEADER + "User interrupted" + bcolors.ENDC
            print(msg)
            with open(self.logfile, "a") as log:
            	log.write(msg)

    def setClustersOutdir(self, clustersOutdir):
        self.clustersOutdir = clustersOutdir
        # Register clusters in each sample 
        for sample in list(self.samples.values()):
            for i in os.listdir(clustersOutdir):
                if os.path.isdir(i) and i == sample.sampleId:
                    sample.registerClustersFromPath(os.path.join(clustersOutdir, sample.sampleId))          

    def velocityAllSamples(self, TE_gtf, geneGTF, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    sampleIndir = os.path.join(self.quantifyOutdir, sample.sampleId)
                    executor.submit(sample.velocity, TE_gtf, geneGTF, sampleIndir)
                executor.shutdown(wait=True)
        except KeyboardInterrupt:
            msg = bcolors.HEADER + "User interrupted" + bcolors.ENDC
            print(msg)
            with open(self.logfile, "a") as log:
            	log.write(msg)

    #def setProcessClustersOutdir(self, processClustersOutdir):
    #    self.processedClustersOutdir = processClustersOutdir

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
            with open(self.logfile, "a") as log:
	            try:
	                # Run job script
	                jobId = jobHandler.runJob("mergeSamples", jobFile, cmd, self.slurm, self.modules)
	                
	                # Print if the job was succesfully submitted
	                msg = jobHandler.sucessSubmit("mergeSamples", self.name, jobId)
	                print(msg)
	                log.write(msg)

	                # Wait for the job to finish and returns an exit code
	                exitCode = jobHandler.waitForJob(jobId)
	            
	                # Print if the job was finished succesfully or not
	                msg = jobHandler.checkExitCodes("mergeSamples", ("Experiment " + self.name), jobId, exitCode)
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
	                            cluster = Cluster(clusterName, os.path.join(outdir, i))
	                            # print(sampleId, cluster)
	                            # print(self.mergeSamples[sampleId].clusters)
	                            self.mergeSamples[sampleId].clusters.append(cluster)
	                    # print([j.clusters for j in self.mergeSamples.values()])
	                return exitCode
	            except:
                    msg = jobHandler.genericError("mergeSamples", self.name)
					print(msg)
					log.write(msg)
	                return
        else:
            subprocess.call(cmd)
            
        self.mergeSamplesOutdir = outdir

    def setMergeSamplesOutdir(self, mergeSamplesOutdir):
    	self.mergeSamplesOutdir = mergeSamplesOutdir
    	for i in os.listdir(mergeSamplesOutdir):
            if(i.endswith(".tsv")):
                clusterName = i.split(".tsv")[0]
                sampleId = clusterName.split("_merged.clusters")[0]
                cluster = Cluster(clusterName, os.path.join(mergeSamplesOutdir, i))
                # print(sampleId, cluster)
                # print(self.mergeSamples[sampleId].clusters)
                self.mergeSamples[sampleId].clusters.append(cluster)
        # print([j.clusters for j in self.mergeSamples.values()])

    def tsvToBamClusters(self, mode, outdir, jobs=1):
        with open(self.logfile, "a") as log:
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
	                return 1

	        except KeyboardInterrupt:
	            msg = bcolors.HEADER + "User interrupted. Finishing tsvToBam for all clusters of all samples before closing." + bcolors.ENDC
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
	                    msg = "Please specify a mode (merged/perSample)"
						print(msg)
						log.write(msg)
	                    return 2

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
	                return 1
	        except KeyboardInterrupt:
	        	msg = bcolors.HEADER + "User interrupted. Finishing filterUMIs for all clusters of all samples before closing." + bcolors.ENDC
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
	                    msg = "Please specify a mode (merged/perSample)"
						print(msg)
						log.write(msg)
	                    return 2

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
	                return 1
	        except KeyboardInterrupt:
	        	msg = bcolors.HEADER + "User interrupted. Finishing bamToFastq for all clusters of all samples before closing." + bcolors.ENDC
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
	                    msg = "Please specify a mode (merged/perSample)"
						print(msg)
						log.write(msg)
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
	                return 1
	        except KeyboardInterrupt:
	        	msg = bcolors.HEADER + "User interrupted. Finishing concatenateLanes for all clusters of all samples before closing." + bcolors.ENDC
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
	                    msg = "Please specify a mode (merged/perSample)"
						print(msg)
						log.write(msg)
	                    return 2

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
	                return 1
	        except KeyboardInterrupt:
	        	msg = bcolors.HEADER + "User interrupted. Finishing mapping for all clusters of all samples before closing." + bcolors.ENDC
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
	                    msg = "Please specify a mode (merged/perSample)"
						print(msg)
						log.write(msg)
	                    return 2

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
	                return 1
	        except KeyboardInterrupt:
	            msg = bcolors.HEADER + "User interrupted. Finishing TEcounts for all clusters of all samples before closing." + bcolors.ENDC
				print(msg)
				log.write(msg)

    def normalizeTECountsCluster(self, mode, jobs=1):
        with open(self.logfile, "a") as log:
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
	                    
	                    msg = jobHandler.sucessSubmit("mergeSamplesNorm", self.name, jobId)
						print(msg)
						log.write(msg)
	                    
	                    exitCode = jobHandler.waitForJob(jobId)
	                    msg = jobHandler.checkExitCodes("mergeSamplesNorm", ("Experiment " + self.name), jobId, exitCode)
						print(msg)
						log.write(msg)

	                else:
	                    subprocess.call(cmd)
	                self.mergeNormalizedTEcountsOutdir = outdir

	                if exitCode:
	                	return 0
		            else:
		                return 1
	            else:
	            	self.normalized_results = []
	                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
	                    for sample in list(self.samples.values()):
	                        sampleIndir = os.path.join(indir, sample.sampleId)
	                        sampleOutdir = os.path.join(outdir, sample.sampleId)
	                        self.normalized_results.append(executor.submit(sample.normalizeTEcounts, sampleIndir, sampleOutdir))
	                        sample.NormalizedOutdir = outdir
	                normalized_exitCodes = [i.result()[1] for i in self.normalized_results]
	            	normalized_allSuccess = all(exitCode == 0 for exitCode in normalized_exitCodes)

	            	if normalized_allSuccess:
	                	return 0
		            else:
		                return 1
	        except KeyboardInterrupt:
	            msg = bcolors.HEADER + "User interrupted" + bcolors.ENDC
				print(msg)
				log.write(msg)

    def processClusters(self, mode, outdir, geneGTF, teGTF, starIndex, jobs=1):
        with open(self.logfile, "a") as log:
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
	            tsvToBam_allSuccess = self.tsvToBamClusters(mode, outdir, jobs)
	            
	            # self.tsvToBam_results = []
	            # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
	            #     for sampleId, sample in samplesDict.items():
	            #         for cluster in sample.clusters:
	            #             bam = os.path.join(self.quantifyOutdir, sampleId, "outs/possorted_genome_bam.bam")
	            #             outdir_sample = os.path.join(outdir, "tsvToBam/", sampleId)
	            #             self.tsvToBam_results.append(executor.submit(cluster.tsvToBam, sampleId, bam, outdir_sample, self.slurm, self.modules))
	            # tsvToBam_exitCodes = [i.result()[1] for i in self.tsvToBam_results]
	            # tsvToBam_allSuccess = all(exitCode == 0 for exitCode in tsvToBam_exitCodes)

	            if tsvToBam_allSuccess:
	                current_instruction = "filterUMIs"
	                filterUMIs_allSuccess = self.filterUMIsClusters(mode, outdir, jobs)

	                # self.filterUMIs_results = []
	                # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
	                #     for sampleId, sample in samplesDict.items():
	                #         for cluster in sample.clusters:
	                #             inbam = os.path.join(outdir, "tsvToBam/", sampleId, (cluster.clusterName + ".bam"))
	                #             outdir_sample = os.path.join(outdir, "filterUMIs/", sampleId)
	                #             #clusters = self.mergeSamplesClusters[sampleId]
	                #             self.filterUMIs_results.append(executor.submit(cluster.filterUMIs, sampleId, inbam, outdir_sample, self.slurm, self.modules))
	                # filterUMIs_exitCodes = [i.result()[1] for i in self.filterUMIs_results]
	                # filterUMIs_allSuccess = all(exitCode == 0 for exitCode in filterUMIs_exitCodes)
	                if filterUMIs_allSuccess:
	                    current_instruction = "bamToFastq"
	                    bamToFastq_allSuccess = self.bamToFastqClusters(mode, outdir, jobs)
	                    
	                    # self.bamToFastq_results = []
	                    # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
	                    #     for sampleId, sample in samplesDict.items():
	                    #         for cluster in sample.clusters:
	                    #             bam = os.path.join(outdir, "filterUMIs/", sampleId, (cluster.clusterName + "_filtered.bam"))
	                    #             outdir_sample = os.path.join(outdir, "bamToFastq/", sampleId)
	                    #             self.bamToFastq_results.append(executor.submit(cluster.bamToFastq, sampleId, bam, outdir_sample, self.slurm, self.modules))
	                    # bamToFastq_exitCodes = [i.result()[1] for i in self.bamToFastq_results]
	                    # bamToFastq_allSuccess = all(exitCode == 0 for exitCode in bamToFastq_exitCodes)
	    
	                    if bamToFastq_allSuccess:
	                        current_instruction = "concatenateLanes"
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
	        	msg = bcolors.HEADER + "User interrupted. Finishing instruction " + current_instruction + " for all clusters of all samples before closing." + bcolors.ENDC
				print(msg)
				log.write(msg)

        