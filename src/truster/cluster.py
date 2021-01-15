#! /bin/python3
import subprocess
import os
from .jobHandler import *
import pysam

class Cluster:
	def __init__(self, clusterName, tsv):
		self.clusterName = clusterName
		self.tsv = tsv
		self.outdirs = {"tsvToBam" : None, "filterUMIs" : None, "bamToFastq" : None, "concatenateLanes" : None, "mapCluster" : None}

	def tsvToBam(self, sampleId, bam, outdir, slurm=None, modules=None):
		try:
			if not os.path.exists("tsvToBam_scripts"):
				os.makedirs("tsvToBam_scripts", exist_ok=True)
			if not os.path.exists(outdir):
				os.makedirs(outdir, exist_ok=True)

			cmd = ["subset-bam", "--bam", bam, "--cell-barcodes", self.tsv, "--out-bam", os.path.join(outdir, (self.clusterName + ".bam"))]

			if slurm != None:
				cmd = ' '.join(cmd)

				jobFile =  os.path.join("tsvToBam_scripts/", (self.clusterName + "_tsvToBam.sh"))
				try:
					jobId = jobHandler.runJob("tsvToBam", jobFile, cmd, slurm, modules)
					print(jobHandler.sucessSubmit("tsvToBam", (" sample " + sampleId + " cluster " +  self.clusterName), jobId))
					exitCode = jobHandler.waitForJob(jobId)
					print(jobHandler.checkExitCodes("tsvToBam", ("Sample " + sampleId + ", cluster " + self.clusterName),jobId, exitCode))
					if exitCode == 0:
						self.outdirs["tsvToBam"] = outdir
					return (jobId, exitCode)
				except:
					print(jobHandler.genericError("tsvToBam", (" sample " + sampleId + " cluster " +  self.clusterName)))
					return
			else:
				subprocess.call(cmd)
		except KeyboardInterrupt:
			print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

	def filterUMIs(self, sampleId, inbam, outdir, slurm=None, modules=None):
		try:
			if not os.path.exists("filterUMIs_scripts"):
				os.makedirs("filterUMIs_scripts", exist_ok=True)
			if not os.path.exists(outdir):
				os.makedirs(outdir, exist_ok=True)
			outbam = os.path.join(outdir, (self.clusterName + "_filtered.bam"))
			cmd = ["python", "filterUMIs.py", "-i", inbam, "-o", outbam]
			if slurm != None:
				cmd = ' '.join(cmd)
				jobFile =  os.path.join("filterUMIs_scripts/", (self.clusterName + "_filterUMIs.sh"))
				try:
					jobId = jobHandler.runJob("filterUMIs", jobFile, cmd, slurm, modules)
					print(jobHandler.sucessSubmit("filterUMIs", (" sample " + sampleId + " cluster " +  self.clusterName), jobId))
					exitCode = jobHandler.waitForJob(jobId)
					print(jobHandler.checkExitCodes("filterUMIs", ("Sample " + sampleId + ", cluster " + self.clusterName),jobId, exitCode))
					if exitCode == 0:
						self.outdirs["filterUMIs"] = outdir
					return (jobId, exitCode)
				except:
					print(jobHandler.genericError("filterUMIs", (" sample " + sampleId + " cluster " +  self.clusterName)))
					return
			else:
				subprocess.call(cmd)
		except KeyboardInterrupt:
			print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

	def bamToFastq(self, sampleId, bam, outdir, slurm=None, modules=None):
		try:
			if not os.path.exists("bamToFastq_scripts"):
				os.makedirs("bamToFastq_scripts", exist_ok=True)
			if not os.path.exists(outdir):
				os.makedirs(outdir, exist_ok=True)
			cmd = ["bamtofastq-1.2.0", bam, (outdir + "/" + self.clusterName)]
			if slurm != None:
				cmd = ' '.join(cmd)
				jobFile =  os.path.join("bamToFastq_scripts/", (self.clusterName + "_bamToFastq.sh"))
				try:
					jobId = jobHandler.runJob("bamToFastq", jobFile, cmd, slurm, modules)
					print(jobHandler.sucessSubmit("bamToFastq", (" sample " + sampleId + " cluster " +  self.clusterName), jobId))
					exitCode = jobHandler.waitForJob(jobId)
					print(jobHandler.checkExitCodes("bamToFastq", ("Sample " + sampleId + ", cluster " + self.clusterName),jobId, exitCode))
					if exitCode == 0:
						self.outdirs["bamToFastq"] = outdir
					return (jobId, exitCode)
				except:
					print(jobHandler.genericError("bamToFastq", (" sample " + sampleId + " cluster " +  self.clusterName)))
					return
			else:
				subprocess.call(cmd)
		except KeyboardInterrupt:
			print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

	def concatenateLanes(self, sampleId, indir, outdir, slurm=None, modules=None):
		try:
			if not os.path.exists("concatenateLanes_scripts"):
				os.makedirs("concatenateLanes_scripts", exist_ok=True)
			if not os.path.exists(outdir):
				os.makedirs(outdir, exist_ok=True)

			cmd = ["cat", os.path.join(indir, self.clusterName, "*/*_R2_001.fastq.gz"), ">", os.path.join(outdir, (self.clusterName + "_R2.fastq.gz"))]
			if slurm != None:
				jobFile =  os.path.join("concatenateLanes_scripts/", (self.clusterName + "_concatenateLanes.sh"))
				try:
					jobId = jobHandler.runJob("concatenateLanes", jobFile, cmd, slurm, modules)
					print(jobHandler.sucessSubmit("concatenateLanes", (" sample " + sampleId + " cluster " +  self.clusterName), jobId))
					exitCode = jobHandler.waitForJob(jobId)
					print(jobHandler.checkExitCodes("concatenateLanes", ("Sample " + sampleId + ", cluster " + self.clusterName),jobId, exitCode))
					if exitCode == 0:
						self.outdirs["concatenateLanes"] = outdir
					return (jobId, exitCode)
				except:
					print(jobHandler.genericError("concatenateLanes", (" sample " + sampleId + " cluster " +  self.clusterName)))
					return
			else:
				subprocess.call([cmd])
		except KeyboardInterrupt:
			print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

	def mapCluster(self, sampleId, fastqdir, outdir, geneGTF, starIndex, RAM, slurm=None, modules=None):
		try:
			if not os.path.exists("mapCluster_scripts"):
				os.makedirs("mapCluster_scripts", exist_ok=True)
			if not os.path.exists(outdir):
				os.makedirs(outdir, exist_ok=True)

			cmd = ["STAR", "--runThreadN", str(slurm["mapCluster"]["tasks-per-node"]), "--readFilesCommand", "gunzip", "-c", "--outSAMattributes", "All", "--outSAMtype", "BAM", "SortedByCoordinate", "--sjdbGTFfile", geneGTF, "--genomeDir", starIndex, "--outFileNamePrefix", (str(os.path.join(outdir, self.clusterName)) + "_") , "--outFilterMultimapNmax", "100", "--limitBAMsortRAM", RAM, "--winAnchorMultimapNmax", "200", "--readFilesIn", os.path.join(fastqdir, (self.clusterName + "_R2.fastq.gz"))]

			if slurm != None:
				cmd = ' '.join(cmd)
				jobFile =  os.path.join("mapCluster_scripts/", (self.clusterName + "_mapCluster.sh"))
				try:
					jobId = jobHandler.runJob("mapCluster", jobFile, cmd, slurm, modules)
					print(jobHandler.sucessSubmit("mapCluster", (" sample " + sampleId + " cluster " +  self.clusterName), jobId))
					exitCode = jobHandler.waitForJob(jobId)
					print(jobHandler.checkExitCodes("mapCluster", ("Sample " + sampleId + ", cluster " + self.clusterName),jobId, exitCode))
					if exitCode == 0:
						self.outdirs["mapCluster"] = outdir
					return (jobId, exitCode)
				except:
					print(jobHandler.genericError("mapCluster", (" sample " + sampleId + " cluster " +  self.clusterName)))
					return
			else:
				subprocess.call(cmd)
		except KeyboardInterrupt:
			print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)

	def TEcount(self, experimentName, sampleId, bam, outdir, geneGTF, teGTF, slurm=None, modules=None):
		try:
			if not os.path.exists("TEcount_scripts"):
				os.makedirs("TEcount_scripts", exist_ok=True)
			if not os.path.exists(outdir):
				os.makedirs(outdir, exist_ok=True)

			cmd = ["TEcount", "-b", bam, "--GTF", geneGTF, "--TE", teGTF, "--format", "BAM", "--stranded", "yes", "--mode", "multi", "--sortByPos", "--project", os.path.join(outdir, (experimentName + "_" + self.clusterName + "_"))]

			if slurm != None:
				cmd = ' '.join(cmd)
				jobFile =  os.path.join("TEcount_scripts/", (self.clusterName + "_TEcount.sh"))
				try:
					jobId = jobHandler.runJob("TEcount", jobFile, cmd, slurm, modules)
					print(jobHandler.sucessSubmit("TEcount", (" sample " + sampleId + " cluster " +  self.clusterName), jobId))
					exitCode = jobHandler.waitForJob(jobId)
					print(jobHandler.checkExitCodes("TEcount", ("Sample " + sampleId + ", cluster " + self.clusterName),jobId, exitCode))
					if exitCode == 0:
						self.outdirs["TEcount"] = outdir
					return (jobId, exitCode)
				except:
					print(jobHandler.genericError("TEcount", (" sample " + sampleId + " cluster " +  self.clusterName)))
					return
			else:
				subprocess.call(cmd)
		except KeyboardInterrupt:
			print(bcolors.HEADER + "User interrupted" + bcolors.ENDC)


