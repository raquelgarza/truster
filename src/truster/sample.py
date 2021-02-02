#! /bin/python3
import subprocess
import os
import json
import sys
from .jobHandler import *
import pysam
import concurrent.futures
from .cluster import Cluster
from .bcolors import Bcolors

class Sample:

    def __init__(self, slurm, modules, logfile, sampleId="", sampleName="", rawPath = ""):
        self.slurm = slurm
        self.modules = modules
        self.sampleId = sampleId
        self.sampleName = sampleName
        self.rawPath = rawPath
        self.logfile = logfile
        # with open(self.logfile, "a") as log:
        #     msg = "Sample " + sampleId + " created\n"
        #     log.write(msg)

    def quantify(self, crIndex, indir, outdir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("counts_scripts/"):
                    os.makedirs("counts_scripts", exist_ok=True)
    
                cmd = ["cellranger count", "--id", self.sampleId, "--transcriptome", crIndex, "--fastqs", indir]
                # If the experiment has a cluster configuration file
                if self.slurm != None:
    
                    cmd = ' '.join([cmd[0], '='.join(cmd[1:3]), '='.join(cmd[3:5]), '='.join(cmd[5:7])])
    
                    jobFile =  os.path.join("counts_scripts/", (self.sampleId + "_counts.sh"))
                    try:
                        # Run a job using runJob from the module which returns a job id
                        jobId = runJob("counts", jobFile, cmd, self.slurm, self.modules)
                        # The latest (or only) "counts" job ran for this sample.
                        # Without this one being completed other functions might not be able to run.
                        msg = sucessSubmit("counts", self.sampleId, jobId)
                        log.write(msg)
                        msg = checkExitCodes("counts", ("Sample " + self.sampleId),jobId, exitCode)
                        log.write(msg)
    
                        exitCode = waitForJob(jobId)
                        
                        if exitCode == 0:
                            subprocess.call("mv", self.sampleId, outdir)
    
                        return exitCode
                    except:
                        # If the job couldnt be created but there is a cluster configuration file...
                        msg = genericError("counts", self.sampleId)
                        log.write(msg)
                        return
                else:
                    # Run locally
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def velocity(self, TE_gtf, geneGTF, indir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                
                cmd = ["velocyto run10x", "-m", TE_gtf, indir, geneGTF]
    
                if self.slurm != None:
    
                    cmd = ' '.join(cmd)
    
                    jobFile =  os.path.join("velocity_scripts/", (self.sampleId + "_velocity.sh"))
                    try:
                        jobId = runJob("velocity", jobFile, cmd, self.slurm, self.modules)
                        msg = sucessSubmit("velocity", self.sampleId, jobId)
                        log.write(msg)

                        exitCode = waitForJob(jobId)
                        msg = checkExitCodes("velocity", ("Sample " + self.sampleId),jobId, exitCode)
                        log.write(msg)
                        
                    except:
                        msg = genericError("velocity", self.sampleId)
                        log.write(msg)
                        return waitForJob(jobId)
                else:
                    subprocess.call(cmd)
                    subprocess.call("mv", self.sampleId, indir)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def plotVelocity(self, loom, indir, outdir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                cwd = os.path.dirname(os.path.realpath(__file__))
                os.path.join(cwd, "py_scripts/plotVelocity.py")
                # print('plotVelocity -l <loom> -n <sample_name> -u <umap> -c <clusters> -o <outdir>')
                cmd = ["python", os.path.join(cwd, "py_scripts/plotVelocity"), "-l", loom, "-n", self.sampleId, "-u", os.path.join(indir, (self.sampleId + "_cell_embeddings.csv")), "-c", os.path.join(indir, (self.sampleId + "_clusters.csv")), "-o", outdir]
    
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

    def emptyClusters(self):
        self.clusters = []
        return

    def getClusters(self, indir, outdir, res = 0.5, percMitochondrial = None, minGenes = None, exclude = None):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("getClusters_scripts"):
                    os.makedirs("getClusters_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                res = str(res)
                cwd = os.path.dirname(os.path.realpath(__file__))
                cmd = [os.path.join(cwd, "r_scripts/get_clusters.R"), "-i", os.path.join(indir, "outs/filtered_feature_bc_matrix"), "-o", outdir, "-s", self.sampleId, "-r", res]

                if exclude != None:
                    cmd.extend(["-e", exclude])
                if percMitochondrial != None:
                    cmd.extend(["-p", str(percMitochondrial)])
                if minGenes != None:
                    cmd.extend(["-m", str(minGenes)])

                if self.slurm != None:
                    cmd = ' '.join(cmd) 
                    jobFile =  os.path.join("getClusters_scripts/", (self.sampleId + "_getClusters.sh"))
                    try:
                        jobId = runJob("getClusters", jobFile, cmd, self.slurm, self.modules)
                        msg = sucessSubmit("getClusters", self.sampleId, jobId)
                        log.write(msg)

                        exitCode = waitForJob(jobId)
                        msg = checkExitCodes("getClusters", ("Sample " + self.sampleId),jobId, exitCode)
                        log.write(msg)
                        
                        if exitCode == 0:
                            self.clusters = [Cluster(j.split(".tsv")[0], os.path.join(outdir, j)) for j in os.listdir(outdir) if j.endswith(".tsv")]
                            self.rdataPath = os.path.join(outdir, (self.sampleId + ".RData"))
                        return exitCode
                    except:
                        msg = genericError("getClusters", self.sampleId)
                        log.write(msg)
                        return
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def registerClustersFromPath(self, path):
        self.clusters = [Cluster(j.split(".tsv")[0], os.path.join(path, j), self.logfile) for j in os.listdir(path) if j.endswith(".tsv")]
        self.rdataPath = os.path.join(path, (self.sampleId + ".RData"))
        msg = "Clusters from " + self.sampleId + " have been registered from " + path + ". Clusters: " + ', '.join([i.clusterName for i in self.clusters]) + "\n"
        return msg

    def normalizeTEcounts(self, indir, outdir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("normalizeTEcounts_scripts"):
                    os.makedirs("normalizeTEcounts_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cwd = os.path.dirname(os.path.realpath(__file__))
                cmd = ["Rscript", os.path.join(cwd, "r_scripts/normalize_TEexpression.R"), "-m", "individual", "-o", outdir, "-i", indir, "-r", self.rdataPath]
                if self.slurm != None:
                    cmd = ' '.join(cmd)
                    jobFile =  os.path.join("normalizeTEcounts_scripts/", (self.sampleId + "_normalizeTEcounts.sh"))
                    try:
                        jobId = runJob("normalizeTEcounts", jobFile, cmd, self.slurm, self.modules)
                        msg = sucessSubmit("normalizeTEcounts", self.sampleId, jobId)
                        log.write(msg)
                        
                        exitCode = waitForJob(jobId)
                        msg = checkExitCodes("normalizeTEcounts", ("Sample " + self.sampleId),jobId, exitCode)
                        log.write(msg)
                        
                        if exitCode == 0:
                            self.normTEcounts = os.path.join(outdir, "TE_normalizedValues_aggregatedByClusters_melted.csv")
                        return exitCode
                    except:
                        msg = genericError("normalizeTEcounts", self.sampleId)
                        log.write(msg)
                        return
                else:
                    subprocess.call(cmd)
                    pass
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)


