#! /bin/python3
import subprocess
from subprocess import PIPE
import pandas as pd
import time
from .bcolors import Bcolors

def runJob(function, job_file, code, slurm, modules):
    with open(job_file, "w") as fout:
        fout.writelines("#!/bin/bash\n")
        try:
            for k,v in slurm[function].items():
                fout.writelines("#SBATCH --" + str(k) + "=" + str(v) + "\n")
        except:
            for k,v in slurm["__default__"].items():
                fout.writelines("#SBATCH --" + str(k) + "=" + str(v) + "\n")
        fout.writelines("module purge\n")
        try:
            for i in modules[function]:
                fout.writelines("module load " + str(i) + "\n")
        except:
            pass
        fout.writelines(code + "\n")
        fout.writelines("module purge\n")

    sbatch_out = subprocess.run([("sbatch " + str(job_file))], shell=True, stdout=PIPE, stderr=PIPE)
    time.sleep(3)
    if sbatch_out.returncode == 0:
        jobId = (sbatch_out.stdout.split()[len(sbatch_out.stdout.split())-1]).decode("utf-8")
        return jobId
    else:
        print(Bcolors.FAIL + "Something went wrong. Please try again" + Bcolors.ENDC)
        return

def cancel(jobId):
    scancel_out = subprocess.run([("scancel " + str(jobId))], shell=True, stdout=PIPE, stderr=PIPE)
    time.sleep(3)
    if scancel_out.returncode == 0:
        print(Bcolors.FAIL + "Job" + str(jobId) + " got cancelled." + Bcolors.ENDC)
        return
    else:
        print(Bcolors.FAIL + "Something went wrong. Please try again" + Bcolors.ENDC)
        return

def checkStatus(jobId):
    job = subprocess.run(("sacct -j " + str(jobId) + " --format=state"), shell=True, stdout=PIPE, stderr=PIPE)
    job = job.stdout.decode("utf-8").split()
    status = job[(len(job)-1)]
    return status.lower()

def onTheWay(jobId):
    status = checkStatus(jobId)
    if status == "running" or status == "pending":
        return True
    else:
        return False

def waitForJob(jobId):
    # print("I'm here!")
    time.sleep(3)
    while onTheWay(jobId):
        time.sleep(3)
    if checkStatus(jobId) == "completed":
        # print(Bcolors.OKGREEN + "Job " + jobId + " has been completed! :-)" + Bcolors.ENDC)
        return 0
    elif checkStatus(jobId) == "cancelled":
        return 1
    elif checkStatus(jobId) == "failed":
        return 2
    else:
        return 3

def checkExitCodes(fun, sampleId_clusterName, jobId, exitCode):
    if exitCode == 0:
        return (Bcolors.OKGREEN + "Job " + jobId + " finished " + fun + " succesfully. " + sampleId_clusterName + Bcolors.ENDC)
    elif exitCode == 1:
        return (Bcolors.FAIL + "Job " + jobId + " was cancelled during " + fun + ". " + sampleId_clusterName + Bcolors.ENDC)
    elif exitCode == 2:
        return (Bcolors.FAIL + "Job " + jobId + " failed during " + fun + ". " + sampleId_clusterName + Bcolors.ENDC)
    elif exitCode == 3:
        return (Bcolors.FAIL + "Something strange happened to job " + jobId + " during " + fun + ". " + sampleId_clusterName + Bcolors.ENDC)

def genericError(fun, info):
    return Bcolors.FAIL + "Something went wrong creating the " + fun + " job for " + info + Bcolors.ENDC

def sucessSubmit(fun, info, jobId):
    return Bcolors.OKBLUE + fun.capitalize() + " for " + info + " has been submitted and has job ID " + jobId + Bcolors.ENDC
