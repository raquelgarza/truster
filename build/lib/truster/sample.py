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

    def __init__(self, slurm, modules, sampleId="", sampleName="", rawPath = ""):
        self.slurm = slurm
        self.modules = modules
        self.sampleId = sampleId
        self.sampleName = sampleName
        self.rawPath = rawPath

    def quantify(self, crIndex, indir, outdir):
        try:
            if not os.path.exists("counts_scripts/"):
                os.makedirs("counts_scripts", exist_ok=True)

            cmd = ["cellranger count", "--id", self.sampleId, "--transcriptome", crIndex, "--fastqs", indir]
            # If the experiment has a cluster configuration file
            if self.slurm != None:

                cmd = ' '.join([cmd[0], '='.join(cmd[1:3]), '='.join(cmd[3:5]), '='.join(cmd[5:7])])

                jobFile =  os.path.join("counts_scripts/", (self.sampleId + "_counts.sh"))
                try:
                    # Run a job using runJob from the jobHandler module which returns a job id
                    jobId = jobHandler.runJob("counts", jobFile, cmd, self.slurm, self.modules)
                    # The latest (or only) "counts" job ran for this sample.
                    # Without this one being completed other functions might not be able to run.
                    print(jobHandler.sucessSubmit("counts", self.sampleId, jobId))

                    exitCode = jobHandler.waitForJob(jobId)
                    print(jobHandler.checkExitCodes("counts", ("Sample " + self.sampleId),jobId, exitCode))
                    if exitCode == 0:
                        subprocess.call("mv", self.sampleId, outdir)

                    return exitCode
                except:
                    # If the job couldnt be created but there is a cluster configuration file...
                    print(jobHandler.genericError("counts", self.sampleId))
                    return
            else:
                # Run locally
                subprocess.call(cmd)
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

    def velocity(self, TE_gtf, geneGTF, indir):
        try:
            if not os.path.exists("velocity_scripts/"):
                os.makedirs("velocity_scripts", exist_ok=True)
            
            cmd = ["velocyto run10x", "-m", TE_gtf, indir, geneGTF]

            if self.slurm != None:

                cmd = ' '.join(cmd)

                jobFile =  os.path.join("velocity_scripts/", (self.sampleId + "_velocity.sh"))
                try:
                    jobId = jobHandler.runJob("velocity", jobFile, cmd, self.slurm, self.modules)
                    print(jobHandler.sucessSubmit("velocity", self.sampleId, jobId))
                    exitCode = jobHandler.waitForJob(jobId)
                    print(jobHandler.checkExitCodes("velocity", ("Sample " + self.sampleId),jobId, exitCode))
                except:
                    print(jobHandler.genericError("velocity", self.sampleId))
                    return jobHandler.waitForJob(jobId)
            else:
                subprocess.call(cmd)
                subprocess.call("mv", self.sampleId, indir)
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

    def emptyClusters(self):
        self.clusters = []
        return

    def getClusters(self, indir, outdir):
        try:
            if not os.path.exists("getClusters_scripts"):
                os.makedirs("getClusters_scripts", exist_ok=True)
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            
            cmd = ["Rscript", "../r_scripts/get_clusters.R", "-i", os.path.join(indir, "outs/filtered_feature_bc_matrix"), "-o", outdir, "-s", self.sampleId]
            if self.slurm != None:
                cmd = ' '.join(cmd) 
                jobFile =  os.path.join("getClusters_scripts/", (self.sampleId + "_getClusters.sh"))
                try:
                    jobId = jobHandler.runJob("getClusters", jobFile, cmd, self.slurm, self.modules)
                    print(jobHandler.sucessSubmit("getClusters", self.sampleId, jobId))
                    exitCode = jobHandler.waitForJob(jobId)
                    print(jobHandler.checkExitCodes("getClusters", ("Sample " + self.sampleId),jobId, exitCode))
                    if exitCode == 0:
                        self.clusters = [clusterTruster(j.split(".tsv")[0], os.path.join(outdir, j)) for j in os.listdir(outdir) if j.endswith(".tsv")]
                        self.rdataPath = os.path.join(outdir, (self.sampleId + ".Rdata"))
                    return exitCode
                except:
                    print(jobHandler.genericError("getClusters", self.sampleId))
                    return
            else:
                subprocess.call(cmd)
        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

    def normalizeTEcounts(self, indir, outdir):
        try:
            if not os.path.exists("normalizeTEcounts_scripts"):
                os.makedirs("normalizeTEcounts_script", exist_ok=True)
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)

                cmd = ["Rscript", "./r_scripts/normalize_TEexpression.R", "-m", "individual", "-o", outdir, "-i", indir, "-r", self.rdataPath]
            if self.slurm != None:
                cmd = ' '.join(cmd)

                jobFile =  os.path.join("normalizeTEcounts_scripts/", (self.sampleId + "_normalizeTEcounts.sh"))
                try:
                    jobId = jobHandler.runJob("normalizeTEcounts", jobFile, cmd, self.slurm, self.modules)
                    print(jobHandler.sucessSubmit("normalizeTEcounts", self.sampleId, jobId))
                    exitCode = jobHandler.waitForJob(jobId)
                    print(jobHandler.checkExitCodes("normalizeTEcounts", ("Sample " + self.sampleId),jobId, exitCode))
                    if exitCode == 0:
                        self.normTEcounts = os.path.join(outdir, "TE_normalizedValues_aggregatedByClusters_melted.csv")
                    return exitCode
                except:
                    print(jobHandler.genericError("normalizeTEcounts", self.sampleId))
                    return
            else:
                subprocess.call(cmd)
                pass

        except KeyboardInterrupt:
            print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)


