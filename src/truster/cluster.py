#! /bin/python3
import subprocess
import os
from .jobHandler import *
import pysam
from .bcolors import Bcolors

class Cluster:
    def __init__(self, cluster_name, tsv, logfile):
        self.cluster_name = cluster_name
        self.tsv = tsv
        self.outdirs = {"tsv_to_bam" : None, "filter_UMIs" : None, "bam_to_fastq" : None, "concatenate_lanes" : None, "map_cluster" : None, "TE_count_unique" : None, "TE_count" : None}
        self.logfile = logfile

    def tsv_to_bam(self, sample_id, bam, outdir, slurm=None, modules=None, dry_run = False):
        if not os.path.exists("tsv_to_bam_scripts"):
            os.makedirs("tsv_to_bam_scripts", exist_ok=True)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)

        with open(self.logfile, "a") as log:
            try:
                cmd = ["subset-bam", "--bam", bam, "--cell-barcodes", self.tsv, "--out-bam", os.path.join(outdir, (self.cluster_name + ".bam"))]
                
                result = run_instruction(cmd = cmd, fun = "tsv_to_bam", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "tsv_to_bam", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["tsv_to_bam"] = outdir
                return result
                    
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def filter_UMIs(self, sample_id, inbam, outdir, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("filter_UMIs_scripts"):
                    os.makedirs("filter_UMIs_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                outbam = os.path.join(outdir, (self.cluster_name + "_filtered.bam"))
                
                cwd = os.path.dirname(os.path.realpath(__file__))
                
                cmd = ["python", os.path.join(cwd, "py_scripts/filterUMIs"), "-i", inbam, "-o", outbam]
                result = run_instruction(cmd = cmd, fun = "filter_UMIs", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "filter_UMIs", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["filter_UMIs"] = outdir
                return result

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def bam_to_fastq(self, sample_id, bam, outdir, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("bam_to_fastq_scripts"):
                    os.makedirs("bam_to_fastq_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                cmd = ["bamtofastq-1.2.0", bam, (outdir + "/" + self.cluster_name)]
                cmd = ["python", os.path.join(cwd, "py_scripts/filterUMIs"), "-i", inbam, "-o", outbam]
                result = run_instruction(cmd = cmd, fun = "bam_to_fastq", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "bam_to_fastq", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["bam_to_fastq"] = outdir
                return result
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def concatenate_lanes(self, sample_id, indir, outdir, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("concatenate_lanes_scripts"):
                    os.makedirs("concatenate_lanes_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cmd = ["cat", os.path.join(indir, self.cluster_name, "*/*_R2_001.fastq.gz"), ">", os.path.join(outdir, (self.cluster_name + "_R2.fastq.gz"))]
                result = run_instruction(cmd = cmd, fun = "concatenate_lanes", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "concatenate_lanes", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["bam_to_fastq"] = outdir
                return result

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def map_cluster(self, sample_id, fastq_dir, outdir, gene_gtf, star_index, RAM, out_tmp_dir = None, unique=False, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("map_cluster_scripts"):
                    os.makedirs("map_cluster_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cmd = ["STAR", "--runThreadN", str(slurm["map_cluster"]["tasks-per-node"]), "--readFilesCommand", "gunzip", "-c", "--outSAMattributes", "All", "--outSAMtype", "BAM", "SortedByCoordinate", "--sjdbGTFfile", str(gene_gtf), "--genomeDir", str(star_index), "--outFileNamePrefix", (str(os.path.join(outdir, self.cluster_name)) + "_") , "--limitBAMsortRAM", str(RAM)]
                if unique:
                    cmd.extend(["--outFilterMultimapNmax", "1", "--outFilterMismatchNoverLmax", "0.03"])
                else:
                    cmd.extend(["--outFilterMultimapNmax", "100", "--winAnchorMultimapNmax", "200"])
                if out_tmp_dir != None:
                    cmd.extend(["--out_tmp_dir", out_tmp_dir])
                cmd.extend(["--readFilesIn", os.path.join(fastq_dir, (self.cluster_name + "_R2.fastq.gz"))])
                result = run_instruction(cmd = cmd, fun = "map_cluster", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "map_cluster", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["map_cluster"] = outdir
                return result

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def TE_count(self, experiment_name, sample_id, bam, outdir, gene_gtf, te_gtf, unique=False, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("TE_count_scripts"):
                    os.makedirs("TE_count_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                
                if unique:
                    cmd = ["featureCounts", "-s", str(1), "-F", "GTF", "-g", "transcript_id", "-a", te_gtf, "-o", os.path.join(outdir, (experiment_name + "_" + self.cluster_name + "_uniqueMap.cntTable")), bam]
                else:
                    cmd = ["TEcount", "-b", bam, "--GTF", gene_gtf, "--TE", te_gtf, "--format", "BAM", "--stranded", "yes", "--mode", "multi", "--sortByPos", "--project", os.path.join(outdir, (experiment_name + "_" + self.cluster_name + "_"))]
                
                if unique:
                    function_name = "TE_count_unique"
                else:
                    function_name = "TE_count"

                result = run_instruction(cmd = cmd, fun = "TE_count", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = function_name, dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs[function_name] = outdir
                return result
                
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)


