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
    def __init__(self, slurm, modules, logfile, sample_id="", sample_name="", raw_path = ""):
        self.slurm = slurm
        self.modules = modules
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.raw_path = raw_path
        self.logfile = logfile

    def quantify(self, cr_index, indir, outdir):
        dry_run = False
        if not os.path.exists("quantify_scripts/"):
            os.makedirs("quantify_scripts", exist_ok=True)

        with open(self.logfile, "a") as log:
            try:
                cmd = ["cellranger count", "--id", self.sample_id, "--transcriptome", cr_index, "--fastqs", indir]
                result = run_instruction(cmd = cmd, fun = "quantify", fun_module = "quantify", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]

                if exit_code == 0:
                    subprocess.call("mv", self.sample_id, outdir)
                    self.quantify_outdir = outdir
                
                return exit_code
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def set_quantification_outdir(self, cellranger_outdir):
        self.quantify_outdir = cellranger_outdir
        with open(self.logfile, "a") as log:
            msg = "Quantification directory for sample" + self.sample_id + " is set to: " + cellranger_outdir + ".\n"
            log.write(msg)

    def velocity(self, te_gtf, gene_gtf, indir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                
                cmd = ["velocyto run10x", "-m", te_gtf, indir, gene_gtf]
                result = run_instruction(cmd = cmd, fun = "velocity", fun_module = "velocity", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]
                
                if exit_code == 0:
                    subprocess.call("mv", self.sample_id, indir)
                
                return exit_code
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def plot_velocity(self, loom, indir, outdir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                cwd = os.path.dirname(os.path.realpath(__file__))
                os.path.join(cwd, "py_scripts/plot_velocity.py")
                # print('plot_velocity -l <loom> -n <sample_name> -u <umap> -c <clusters> -o <outdir>')
                cmd = ["python", os.path.join(cwd, "py_scripts/plot_velocity"), "-l", loom, "-n", self.sample_id, "-u", os.path.join(indir, (self.sample_id + "_cell_embeddings.csv")), "-c", os.path.join(indir, (self.sample_id + "_clusters.csv")), "-o", outdir]
                result = run_instruction(cmd = cmd, fun = "plot_velocity", fun_module = "plot_velocity", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]
                return exit_code
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def empty_clusters(self):
        self.clusters = []
        return

    def get_clusters(self, outdir, group_name, res = 0.5, perc_mitochondrial = None, min_genes = None, max_genes = 7000, normalization_method = "LogNormalize", max_size=500):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("get_clusters_scripts"):
                    os.makedirs("get_clusters_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                res = str(res)
                max_size = str(max_size)
                max_genes = str(max_genes)
                
                cwd = os.path.dirname(os.path.realpath(__file__))
                cmd = [os.path.join(cwd, "r_scripts/get_clusters.R"), "-i", os.path.join(self.quantify_outdir, "outs/filtered_feature_bc_matrix"), "-o", outdir, "-s", self.sample_id, "-r", res, "-n", normalization_method, "-S", max_size, "-M", max_genes]

                if perc_mitochondrial != None:
                    cmd.extend(["-p", str(perc_mitochondrial)])
                if min_genes != None:
                    cmd.extend(["-m", str(min_genes)])

                result = run_instruction(cmd = cmd, fun = "get_clusters", fun_module = "get_clusters", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.clusters = [Cluster(j.split(".tsv")[0], os.path.join(outdir, j), self.logfile) for j in os.listdir(outdir) if j.endswith(".tsv")]
                    self.rdata_path = os.path.join(outdir, (self.sample_id + ".rds"))
                return exit_code

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def register_clusters_from_path(self, path):
        self.clusters = [Cluster(j.split(".tsv")[0], os.path.join(path, j), self.logfile) for j in os.listdir(path) if j.endswith(".tsv")]
        self.rdata_path = os.path.join(path, (self.sample_id + ".rds"))
        msg = "Clusters from " + self.sample_id + " have been registered from " + path + ". Clusters: " + ', '.join([i.cluster_name for i in self.clusters]) + "\n"
        return msg

    def normalize_TE_counts(self, indir, outdir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("normalize_TE_counts_scripts"):
                    os.makedirs("normalize_TE_counts_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cwd = os.path.dirname(os.path.realpath(__file__))
                cmd = ["Rscript", os.path.join(cwd, "r_scripts/normalize_TEexpression.R"), "-m", "individual", "-o", outdir, "-i", indir, "-r", self.rdata_path, "-n", self.sample_id]
                
                result = run_instruction(cmd = cmd, fun = "normalize_TE_expression", fun_module = "normalize_TE_expression", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]
                        
                if exit_code == 0:
                    self.norm_TE_counts = os.path.join(outdir, "TE_normalizedValues_aggregatedByClusters_melted.csv")
                return exit_code

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)


