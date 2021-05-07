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
        self.outdirs = {"tsv_to_bam" : None, "filter_UMIs" : None, "bam_to_fastq" : None, "concatenate_lanes" : None, "map_cluster" : None}
        self.logfile = logfile

    def tsv_to_bam(self, sample_id, bam, outdir, slurm=None, modules=None):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("tsv_to_bam_scripts"):
                    os.makedirs("tsv_to_bam_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cmd = ["subset-bam", "--bam", bam, "--cell-barcodes", self.tsv, "--out-bam", os.path.join(outdir, (self.cluster_name + ".bam"))]
    
                if slurm != None:
                    cmd = ' '.join(cmd)
    
                    job_file =  os.path.join("tsv_to_bam_scripts/", (self.cluster_name + "_tsv_to_bam.sh"))
                    try:
                        job_id = run_job("tsv_to_bam", job_file, cmd, slurm, modules)
                        msg = sucess_submit("tsv_to_bam", (" sample " + sample_id + " cluster " +  self.cluster_name), job_id)
                        log.write(msg)

                        exit_code = wait_for_job(job_id)
                        msg = check_exit_codes("tsv_to_bam", ("Sample " + sample_id + ", cluster " + self.cluster_name),job_id, exit_code)
                        log.write(msg)
                        if exit_code == 0:
                            self.outdirs["tsv_to_bam"] = outdir
                        return (job_id, exit_code)
                    except:
                        msg = generic_error("tsv_to_bam", (" sample " + sample_id + " cluster " +  self.cluster_name))
                        log.write(msg)
                        return
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def filter_UMIs(self, sample_id, inbam, outdir, slurm=None, modules=None):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("filter_UMIs_scripts"):
                    os.makedirs("filter_UMIs_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                outbam = os.path.join(outdir, (self.cluster_name + "_filtered.bam"))
                
                cwd = os.path.dirname(os.path.realpath(__file__))
                
                cmd = ["python", os.path.join(cwd, "py_scripts/filterUMIs"), "-i", inbam, "-o", outbam]
                if slurm != None:
                    cmd = ' '.join(cmd)
                    job_file =  os.path.join("filter_UMIs_scripts/", (self.cluster_name + "_filter_UMIs.sh"))
                    try:
                        job_id = run_job("filter_UMIs", job_file, cmd, slurm, modules)
                        msg = sucess_submit("filter_UMIs", (" sample " + sample_id + " cluster " +  self.cluster_name), job_id)
                        log.write(msg)

                        exit_code = wait_for_job(job_id)
                        msg = check_exit_codes("filter_UMIs", ("Sample " + sample_id + ", cluster " + self.cluster_name),job_id, exit_code)
                        log.write(msg)
                        if exit_code == 0:
                            self.outdirs["filter_UMIs"] = outdir
                        return (job_id, exit_code)
                    except:
                        msg = generic_error("filter_UMIs", (" sample " + sample_id + " cluster " +  self.cluster_name))
                        log.write(msg)
                        return
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def bam_to_fastq(self, sample_id, bam, outdir, slurm=None, modules=None):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("bam_to_fastq_scripts"):
                    os.makedirs("bam_to_fastq_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                cmd = ["bamtofastq-1.2.0", bam, (outdir + "/" + self.cluster_name)]
                if slurm != None:
                    cmd = ' '.join(cmd)
                    job_file =  os.path.join("bam_to_fastq_scripts/", (self.cluster_name + "_bam_to_fastq.sh"))
                    try:
                        job_id = run_job("bam_to_fastq", job_file, cmd, slurm, modules)
                        msg = sucess_submit("bam_to_fastq", (" sample " + sample_id + " cluster " +  self.cluster_name), job_id)
                        log.write(msg)

                        exit_code = wait_for_job(job_id)
                        msg = check_exit_codes("bam_to_fastq", ("Sample " + sample_id + ", cluster " + self.cluster_name),job_id, exit_code)
                        log.write(msg)
                        
                        if exit_code == 0:
                            self.outdirs["bam_to_fastq"] = outdir
                        return (job_id, exit_code)
                    except:
                        msg = generic_error("bam_to_fastq", (" sample " + sample_id + " cluster " +  self.cluster_name))
                        log.write(msg)
                        return
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def concatenate_lanes(self, sample_id, indir, outdir, slurm=None, modules=None):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("concatenate_lanes_scripts"):
                    os.makedirs("concatenate_lanes_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cmd = ["cat", os.path.join(indir, self.cluster_name, "*/*_R2_001.fastq.gz"), ">", os.path.join(outdir, (self.cluster_name + "_R2.fastq.gz"))]
                if slurm != None:
                    cmd = ' '.join(cmd)
                    job_file =  os.path.join("concatenate_lanes_scripts/", (self.cluster_name + "_concatenate_lanes.sh"))
                    try:
                        job_id = run_job("concatenate_lanes", job_file, cmd, slurm, modules)
                        msg = sucess_submit("concatenate_lanes", (" sample " + sample_id + " cluster " +  self.cluster_name), job_id)
                        log.write(msg)

                        exit_code = wait_for_job(job_id)
                        msg = check_exit_codes("concatenate_lanes", ("Sample " + sample_id + ", cluster " + self.cluster_name),job_id, exit_code)
                        log.write(msg)
                        if exit_code == 0:
                            self.outdirs["concatenate_lanes"] = outdir
                        return (job_id, exit_code)
                    except:
                        msg = generic_error("concatenate_lanes", (" sample " + sample_id + " cluster " +  self.cluster_name))
                        log.write(msg)
                        return
                else:
                    subprocess.call([cmd])
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def map_cluster(self, sample_id, fastq_dir, outdir, gene_gtf, star_index, RAM, out_tmp_dir = None, unique=False, slurm=None, modules=None):
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
                if slurm != None:
                    cmd = ' '.join(cmd)
                    job_file =  os.path.join("map_cluster_scripts/", (self.cluster_name + "_map_cluster.sh"))
                    try:
                        job_id = run_job("map_cluster", job_file, cmd, slurm, modules)
                        msg = sucess_submit("map_cluster", (" sample " + sample_id + " cluster " +  self.cluster_name), job_id)
                        log.write(msg)

                        exit_code = wait_for_job(job_id)
                        msg = check_exit_codes("map_cluster", ("Sample " + sample_id + ", cluster " + self.cluster_name),job_id, exit_code)
                        log.write(msg)
                        
                        if exit_code == 0:
                            self.outdirs["map_cluster"] = outdir
                        return (job_id, exit_code)
                    except:
                        msg = generic_error("map_cluster", (" sample " + sample_id + " cluster " +  self.cluster_name))
                        log.write(msg)
                        return
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def TE_count(self, experiment_name, sample_id, bam, outdir, gene_gtf, te_gtf, unique=False, slurm=None, modules=None):
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
    
                if slurm != None:
                    cmd = ' '.join(cmd)
                    job_file =  os.path.join("TE_count_scripts/", (self.cluster_name + "_TE_count.sh"))
                    try:
                        if unique:
                            function_name = "TE_count_unique"
                        else:
                            function_name = "TE_count"
                        job_id = run_job(function_name, job_file, cmd, slurm, modules)
                        msg = sucess_submit(function_name, (" sample " + sample_id + " cluster " +  self.cluster_name), job_id)
                        log.write(msg)
                        
                        exit_code = wait_for_job(job_id)
                        msg = check_exit_codes(function_name, ("Sample " + sample_id + ", cluster " + self.cluster_name),job_id, exit_code)
                        log.write(msg)
                        
                        if exit_code == 0:
                            self.outdirs["TE_count"] = outdir
                        return (job_id, exit_code)
                    except:
                        msg = generic_error(function_name, (" sample " + sample_id + " cluster " +  self.cluster_name))
                        log.write(msg)
                        return
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)


