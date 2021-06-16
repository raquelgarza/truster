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

    def __init__(self, name="", slurm_path=None, modules_path=None):
        self.name = name
        self.slurm_path = slurm_path
        self.modules_path = modules_path
        self.samples = {}
        self.logfile = self.name + ".log"

        with open(self.logfile, "w+") as log:
            msg = "Project " + self.name + " created.\n"
            log.write(msg)
            if slurm_path != None:
                try:
                    with open(slurm_path, "r") as config_file:
                        self.slurm = json.load(config_file)
                        msg = "Configuration file loaded.\n"
                        log.write(msg)
                except FileNotFoundError:
                    msg = Bcolors.FAIL + "Error: Cluster configuration file not found" + Bcolors.ENDC + "\n"
                    print(msg)
                    log.write(msg)
                    return 1

            if modules_path != None:
                try:
                    with open(modules_path, "r") as modules_file:
                        self.modules = json.load(modules_file)
                        msg = "Software modules json loaded.\n"
                        log.write(msg)
                except FileNotFoundError:
                    msg = Bcolors.FAIL + "Error: Module configuration file not found" + Bcolors.ENDC + "\n"
                    print(msg)
                    log.write(msg)
                    return 1

    def register_sample(self, sample_id = "", sample_name = "", raw_path = ""):
        new_sample = {sample_id : Sample(slurm=self.slurm, modules = self.modules, sample_id = sample_id, sample_name = sample_name, raw_path = raw_path, logfile = self.logfile)}
        self.samples = {**self.samples, **new_sample}
        
        with open(self.logfile, "a") as log:
            msg = "Sample " + sample_id + " registered.\n"
            log.write(msg)
        # self.samples.extend(sampleTruster(slurm=self.slurm, modules = self.modules, sample_id = sample_id, sample_name = sample_name))

    def unregister_sample(self, sample_id):
        self.samples.pop(sample_id)
        with open(self.logfile, "a") as log:
            msg = "Sample " + sample_id + " unregistered.\n"
            log.write(msg)

    def register_samples_from_path(self, indir, folder_names_as_sample_ids=True):
        def path_to_samples(path):
            # Returns a list of lists containing the name of the last directory before the lane files and a file
            lane_files = [(os.path.join(directories_fullpath, file)).split("/")[-2:] for directories_fullpath, directories_names, file_names in os.walk(os.path.expanduser(path)) for file in file_names]
            # The name of the file will be the sample name and the name of the directory the sample id
            sample_id_name = {lane[0] : lane[1].split("_L00")[0] for lane in lane_files}
            new_samples = {sample_id : Sample(slurm=self.slurm, modules = self.modules, sample_id = sample_id, sample_name = sample_name, raw_path = os.path.join(path, sample_id), logfile = self.logfile) for sample_id, sample_name in sample_id_name.items()}
            return {**self.samples, **new_samples}

        with open(self.logfile, "a") as log:
            if(folder_names_as_sample_ids):
                if isinstance(indir, list):
                    for path in indir:
                        if os.path.isdir(path):
                            self.samples = path_to_samples(path)
                            msg = "Registered from " + path + "\n"
                            log.write(msg)
                elif isinstance(indir, str):
                    if os.path.isdir(indir):
                        self.samples = path_to_samples(indir)
                        msg = "Registered from " + path + "\n"
                        log.write(msg)
                else:
                    msg = "When registering samples from path: indir does not exist (Check " + path + ")\n"
                    log.write(msg)
                with open(self.logfile, "a") as log:
                    msg = "Registered samples: " + str(', '.join([sample.sample_id for sample in list(self.samples.values())]) + ".\n")
                    log.write(msg)

    def quantify(self, cr_index, outdir, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in self.samples.values():
                    sample_indir = sample.raw_path
                    sample_outdir = os.path.join(outdir, sample.sample_id)
                    executor.submit(sample.quantify, cr_index, sample_indir, sample_outdir)
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + ".\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    def set_quantification_outdir(self, sample_id, cellranger_outdir):
        self.samples[sample_id].set_quantification_outdir(cellranger_outdir)

    def get_clusters_all_samples(self, outdir, groups = None, res = 0.5, perc_mitochondrial = None, min_genes = None, max_genes = None, normalization_method = "LogNormalize", max_size = 500, jobs=1):
        with open(self.logfile, "a") as log:
            try:
                # if groups != None:
                #     if isinstance(groups, list) and all([isinstance(i, list) for i in groups]):

                #     group = list(self.merge_samples.keys())
                #     merge_cluster_per_group_out.append(merge_cluster_per_group(group, group_names))
                # else:
                #     if len(groups) == len(group_names):
                #         for i in range(0,len(groups)):
                #             group = groups[i]
                #             group_name = group_names[i]
                #             log.write("Running merge_clusters for samples " + str(group) + ", members of group " + group_name)
                #             if all([sample_id in self.merge_samples.keys() for sample_id in group]):
                #                 log.write("All samples of group " + group_name + " have been found registered.")
                #                 merge_cluster_per_group_out.append(merge_cluster_per_group(group, group_name))
                #             else:
                #                 msg = "Not all the sample ids were found in self.merge_samples. Are you sure you are passing sample ids and not sample names?"
                #                 log.write(msg)
                #     else:
                #         msg = "The list of groups and the list of the name of the groups should be of equal length."
                #         log.write(msg)

                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample in list(self.samples.values()):
                        # if os.path.isdir(sample.quantify_outdir):
                        #     sample_indir = sample.quantify_outdir
                        # else:
                        #     msg = "Error: File not found. Please make sure that " + sample.quantify_outdir + " exists.\n"
                        #     log.write(msg)
                        #     return 1
                        sample_outdir = os.path.join(outdir, sample.sample_id)
                        res = str(res)
                        max_size = str(max_size)
                        min_genes = str(min_genes)
                        max_genes = str(max_genes)
                        
                        msg = "Clustering " + sample.sample_id + " using all cells.\n"
                        executor.submit(sample.get_clusters, sample_outdir, group_name, res, perc_mitochondrial, min_genes, max_genes, normalization_method, max_size)
                        self.clusters_outdir = outdir
                        log.write(msg)
                
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + ".\n"
                log.write(msg)

    def set_clusters_outdir(self, clusters_outdir):
        with open(self.logfile, "a") as log:
            self.clusters_outdir = clusters_outdir
            # Register clusters in each sample 
            for sample in list(self.samples.values()):
                for i in os.listdir(clusters_outdir):
                    if os.path.isdir(os.path.join(clusters_outdir, i)) and  i == sample.sample_id:
                        msg = sample.register_clusters_from_path(os.path.join(clusters_outdir, sample.sample_id)) 
                        log.write(msg)
                    elif os.path.isdir(os.path.join(clusters_outdir, i)) and  i == sample.sample_name:
                        msg = sample.register_clusters_from_path(os.path.join(clusters_outdir, sample.sample_id)) 
                        log.write(msg)
            msg = "The directory for clusters of individual samples is set to: " + clusters_outdir + ".\n"
            log.write(msg)

    def velocity_all_samples(self, te_gtf, gene_gtf, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    if os.path.isdir(sample.quantify_outdir):
                            sample_indir = sample.quantify_outdir
                    else:
                        msg = "Error: File not found. Please make sure that " + sample.quantify_outdir + " exists.\n"
                        log.write(msg)
                        return 1
                    executor.submit(sample.velocity, te_gtf, gene_gtf, sample_indir)
                executor.shutdown(wait=True)
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + "\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    def plot_velocity_all_samples(self, jobs=1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for sample in list(self.samples.values()):
                    if os.path.isdir(sample.quantify_outdir):
                            sample_indir = sample.quantify_outdir
                    else:
                        msg = "Error: File not found. Please make sure that " + sample.quantify_outdir + " exists.\n"
                        log.write(msg)
                        return 1
                    loom = os.path.join(sample.quantify_outdir, "velocyto", (sample.sample_id + ".loom"))
                    sample_outdir = os.path.join(sample.quantify_outdir, "velocyto", "plots")
                    executor.submit(sample.plot_velocity, loom, sample_indir, sample_outdir)
                executor.shutdown(wait=True)
        except KeyboardInterrupt:
            msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n" + "\n"
            with open(self.logfile, "a") as log:
                log.write(msg)

    def plot_velocity_merged(self, loom, names):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                cwd = os.path.dirname(os.path.realpath(__file__))
                os.path.join(cwd, "py_scripts/plot_velocity.py")
                outdir = self.outdir_merged_clusters
                # print('plot_velocity -l <loom> -n <sample_name> -u <umap> -c <clusters> -o <outdir>')
                cmd = ["python", os.path.join(cwd, "py_scripts/plot_velocity"), "-l", ','.join(loom), "-n", ','.join(names), "-u", ','.join([os.path.join(outdir, (name + "_cell_embeddings.csv")) for name in names]), "-c", ','.join([os.path.join(outdir, (name + "_clusters.csv")) for name in names]), "-o", outdir]
    
                if self.slurm != None:
                    cmd = ' '.join(cmd)
                    job_file =  os.path.join("velocity_scripts/", (self.sample_id + "_plot_velocity.sh"))
                    try:
                        job_id = run_job("plot_velocity", job_file, cmd, self.slurm, self.modules)
                        msg = sucess_submit("plot_velocity", self.sample_id, job_id)
                        log.write(msg)

                        exit_code = wait_for_job(job_id)
                        msg = check_exit_codes("plot_velocity", ("Sample " + self.sample_id), job_id, exit_code)
                        log.write(msg)
                        
                    except:
                        msg = generic_error("plot_velocity", self.sample_id)
                        log.write(msg)
                        return wait_for_job(job_id)
                else:
                    subprocess.call(cmd)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)
    
    def merge_samples(self, outdir, normalization_method, integrate_samples = "FALSE", max_size=500):
        # Rscript {input.script} -i {rdata} -n {samplenames} -o {params.outpath}
        # Paths to RData files
        with open(self.logfile, "a") as log:
            samples_seurat_rds = [sample.rdata_path for sample in list(self.samples.values())]
            
            # Sample ids
            samples_ids = [sample.sample_id for sample in list(self.samples.values())]
            
            msg = "Merging samples to produce a combined clustering.\n"
            log.write(msg)
            # If we haven't made the merge before, create a directory to store the scripts needed to do so
            if not os.path.exists("merge_samples_scripts"):
                os.makedirs("merge_samples_scripts", exist_ok=True)
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            
            max_size = str(max_size)
            # Run script ../r_scripts/merge_samples.R with input (-i) of the RData paths
            # and output (-o) of the output directory desired, -s for sample ids,
            # and -e for sample names used in cellranger
            cwd = os.path.dirname(os.path.realpath(__file__))
            cmd = ["Rscript", os.path.join(cwd, "r_scripts/merge_samples.R"), "-i", ','.join(samples_seurat_rds), "-o", outdir, "-s", ','.join(samples_ids), "-e", self.name, "-n", normalization_method, "-S", max_size, "-I", integrate_samples]
            
            # If we are on a server with slurm, use the configuration file to send the jobs 
            if self.slurm_path != None:    
                cmd = ' '.join(cmd)
                
                # This command will be stored along with the slurm configurations at
                # merge_samples_scripts/sample_merge_samples.sh
                job_file =  "merge_samples_scripts/" + self.name + "_merge_samples.sh"
                try:
                    # Run job script
                    job_id = run_job("merge_samples", job_file, cmd, self.slurm, self.modules)
                    
                    # Print if the job was succesfully submitted
                    msg = sucess_submit("merge_samples", self.name, job_id)
                    print(msg)
                    log.write(msg)
                    # Wait for the job to finish and returns an exit code
                    exit_code = wait_for_job(job_id)
                    print("Exit code: " + str(exit_code))
                    # Print if the job was finished succesfully or not
                    msg = check_exit_codes("merge_samples", ("Experiment " + self.name), job_id, exit_code)
                    print(msg)
                    log.write(msg)
                    
                    # If it finished succesfully then 
                    if exit_code == 0:
                
                        # For each of the registered samples
                        self.merge_samples = copy.deepcopy(self.samples)
    
                        # Empty the clusters bc we made new ones (shared/merged)
                        for k,v in self.merge_samples.items():
                            v.empty_clusters()
                    
                        msg = "Emptied clusters"
                        print(msg)
                        log.write(msg)
                    
                        # print([j.clusters for j in self.merge_samples.values()])
    
                        # Make a dictionary of the same sort as the registered samples
                        # for example {sample1 : [cluster1, cluster2]}
                        # with the clusters that we created in the outdir we passed to R
                        # They all have the words "merged.clusters", after that is the number
                        # Before that is the sample id, which we can use as a key in the 
                        # merge_samples_clusters dictionary and just append the cluster objects
                        # To the empty list we now have
                        for i in os.listdir(outdir):
                            if(i.endswith(".tsv")):
                                cluster_name = i.split(".tsv")[0]
                                sample_id = cluster_name.split("_merged.clusters")[0]
                                cluster = Cluster(cluster_name = cluster_name, tsv = os.path.join(outdir, i), logfile = self.logfile)
                                self.merge_samples[sample_id].clusters.append(cluster)
                        # print([j.clusters for j in self.merge_samples.values()])
                        self.merge_samples_outdir = outdir
                    return exit_code
                except:
                    msg = generic_error("merge_samples", self.name)
                    print(msg)
                    log.write(msg)
                    return
            else:
                subprocess.call(cmd)
                self.merge_samples_outdir = outdir

    def set_merge_samples_outdir(self, merge_samples_outdir):
        self.merge_samples_outdir = merge_samples_outdir
        # For each of the registered samples
        self.merge_samples = copy.deepcopy(self.samples)
    
        # Empty the clusters bc we made new ones (shared/merged)
        for k,v in self.merge_samples.items():
            v.empty_clusters()

        for i in os.listdir(merge_samples_outdir):
            if(i.endswith(".tsv")):
                cluster_name = i.split(".tsv")[0]
                sample_id = cluster_name.split("_merged.clusters")[0]
                cluster = Cluster(cluster_name = cluster_name, tsv = os.path.join(merge_samples_outdir, i), logfile = self.logfile)

                self.merge_samples[sample_id].clusters.append(cluster)

        with open(self.logfile, "a") as log:
            msg = "The directory for clusters of combined samples is set to: " + merge_samples_outdir + ".\n\n"
            log.write(msg)

            registered_clusters = [sample.clusters for sample_id,sample in self.merge_samples.items()]
            names_registered_clusters = ', '.join([cluster.cluster_name for list_of_clusters in registered_clusters for cluster in list_of_clusters])
            msg = "Registered merge_samples clusters per sample: " + names_registered_clusters + "\n\n"
            log.write(msg)

    def tsv_to_bam_clusters(self, mode, outdir, jobs=1):
        print("Running tsv_to_bam with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samples_dict = self.merge_samples
                else:
                    if mode == "per_sample":
                        samples_dict = self.samples
                    else:
                        msg = "Please specify a mode (merged/per_sample).\n"
                        print(msg)
                        log.write(msg)
                        return 2
                    
                self.tsv_to_bam_results = []
                msg = "Extracting cell barcodes from BAM files.\n"
                log.write(msg)

                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample_id, sample in samples_dict.items():
                        for cluster in sample.clusters:
                            if os.path.isdir(sample.quantify_outdir):
                                bam = os.path.join(sample.quantify_outdir, "outs/possorted_genome_bam.bam")
                            else:
                                msg = "Error: File not found. Please make sure that " + sample.quantify_outdir + " exists.\n"
                                log.write(msg)
                                return 1
                            outdir_sample = os.path.join(outdir, "tsv_to_bam/", sample_id)
                            self.tsv_to_bam_results.append(executor.submit(cluster.tsv_to_bam, sample_id, bam, outdir_sample, self.slurm, self.modules))

                tsv_to_bam_exit_codes = [i.result()[1] for i in self.tsv_to_bam_results]
                tsv_to_bam_all_success = all(exit_code == 0 for exit_code in tsv_to_bam_exit_codes)
                
                if tsv_to_bam_all_success:
                    msg = "\ntsv_to_bam finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\ntsv_to_bam did not finished succesfully for all samples\n"
                    log.write(msg)
                    return False

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing tsv_to_bam for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def filter_UMIs_clusters(self, mode, outdir, jobs=1):
        print("Running filter_UMIs with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samples_dict = self.merge_samples
                else:
                    if mode == "per_sample":
                        samples_dict = self.samples
                    else:
                        msg = "Please specify a mode (merged/per_sample).\n"
                        print(msg)
                        log.write(msg)
                        return 2

                self.filter_UMIs_results = []
                msg = "Extracting cell barcodes from BAM files.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample_id, sample in samples_dict.items():
                        for cluster in sample.clusters:
                            inbam = os.path.join(outdir, "tsv_to_bam/", sample_id, (cluster.cluster_name + ".bam"))
                            outdir_sample = os.path.join(outdir, "filter_UMIs/", sample_id)
                            self.filter_UMIs_results.append(executor.submit(cluster.filter_UMIs, sample_id, inbam, outdir_sample, self.slurm, self.modules))
                            
                filter_UMIs_exit_codes = [i.result()[1] for i in self.filter_UMIs_results]
                filter_UMIs_all_success = all(exit_code == 0 for exit_code in filter_UMIs_exit_codes)
                
                if filter_UMIs_all_success:
                    msg = "\nfilter_UMIs finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nfilter_UMIs did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing filter_UMIs for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)
 
    def bam_to_fastq_clusters(self, mode, outdir, jobs=1):
        print("Running bam_to_fastq with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samples_dict = self.merge_samples
                else:
                    if mode == "per_sample":
                        samples_dict = self.samples
                    else:
                        msg = "Please specify a mode (merged/per_sample).\n"
                        print(msg)
                        log.write(msg)
                        return 2

                self.bam_to_fastq_results = []
                msg = "Converting BAM to FastQ files.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample_id, sample in samples_dict.items():
                        for cluster in sample.clusters:
                            bam = os.path.join(outdir, "filter_UMIs/", sample_id, (cluster.cluster_name + "_filtered.bam"))
                            outdir_sample = os.path.join(outdir, "bam_to_fastq/", sample_id)
                            self.bam_to_fastq_results.append(executor.submit(cluster.bam_to_fastq, sample_id, bam, outdir_sample, self.slurm, self.modules))
                            
                bam_to_fastq_exit_codes = [i.result()[1] for i in self.bam_to_fastq_results]
                bam_to_fastq_all_success = all(exit_code == 0 for exit_code in bam_to_fastq_exit_codes)

                if bam_to_fastq_all_success:
                    msg = "\nbamToFastq finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nbamToFastq did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing bam_to_fastq for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def concatenate_lanes_clusters(self, mode, outdir, jobs=1):
        print("Running concatenate_lanes with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                if mode == "merged":
                    samples_dict = self.merge_samples
                else:
                    if mode == "per_sample":
                        samples_dict = self.samples
                    else:
                        msg = "Please specify a mode (merged/per_sample).\n"
                        print(msg)
                        log.write(msg)
                        return

                self.concatenate_lanes_results = []
                msg = "Concatenating FastQ files.\n"
                log.write(msg)
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    for sample_id, sample in samples_dict.items():
                        for cluster in sample.clusters:
                            indir = os.path.join(outdir, "bam_to_fastq/", sample_id)
                            outdir_sample = os.path.join(outdir, "concatenate_lanes/", sample_id)
                            self.concatenate_lanes_results.append(executor.submit(cluster.concatenate_lanes, sample_id, indir, outdir_sample, self.slurm, self.modules))
                            
                concatenate_lanes_exit_codes = [i.result()[1] for i in self.concatenate_lanes_results]
                concatenate_lanes_all_success = all(exit_code == 0 for exit_code in concatenate_lanes_exit_codes)

                if concatenate_lanes_all_success:
                    msg = "\nconcatenate_lanes finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nconcatenate_lanes did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing concatenate_lanes for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def merge_clusters(self, outdir, groups = {"merged_cluster": "all"}):
        with open(self.logfile, "a") as log:
            def merge_cluster_per_group(group, group_name):
                merge_samples_lists_clusters = [samples_clusters.clusters for sample_id, samples_clusters in self.merge_samples.items() if sample_id in group]
                merge_samples_clusters = [cluster.cluster_name.split("merged.clusters_")[1] for merge_samples_list_clusters in merge_samples_lists_clusters for cluster in merge_samples_list_clusters]
                
                tsvs = { key : set() for key in merge_samples_clusters }
                for sample_id, samples_clusters in self.merge_samples.items():
                    if sample_id in group:
                        for cluster in samples_clusters.clusters:
                            cluster_num = cluster.cluster_name.split("merged.clusters_")[1]
                            tsvs[cluster_num].add(cluster.tsv)

                merge_clusters = dict.fromkeys(tsvs.keys())
                for cluster_num, tsv in tsvs.items():
                    merge_clusters[cluster_num] = Cluster(cluster_name = (group_name + "_" + str(cluster_num)), tsv = list(tsvs[cluster_num]), logfile = self.logfile)
                
                self.merged_clusters_results = []
                # with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                for cluster_num in merge_clusters.keys():
                    outdir_concatenat_lanes = os.path.join(outdir, "concatenate_lanes")
                    outfile = os.path.join(outdir_merged_clusters, (group_name + "_" + str(cluster_num) + "_R2.fastq.gz"))
                    
                    cluster_fastqs = []
                    for dirpath, subdirs, files in os.walk(outdir_concatenat_lanes):
                        for file in files:
                            if file.split("_merged.clusters_")[1].split("_R2.fastq.gz")[0] == cluster_num and file.split("_merged.clusters_")[0] in group:
                                cluster_fastqs.append(os.path.join(dirpath, file))

                    cmd = cluster_fastqs
                    cmd.insert(0, "cat")
                    
                    log.write("Running " + ' '.join(cmd) + "\n\n\n")
                    
                    with open(outfile, "w") as fout:
                        self.merged_clusters_results.append(subprocess.call(cmd, shell = False, stdout=fout, universal_newlines=True))
                
                self.merge_samples_groups[group_name] = merge_clusters
                self.outdir_merged_clusters_groups[group_name] = outdir_merged_clusters
                log.write("self.merge_samples_groups[" + group_name + "] is now " + str(merge_clusters) + "\n\n\n")
                
                merged_clusters_all_success = all(exit_code == 0 for exit_code in self.merged_clusters_results)

                if merged_clusters_all_success:
                    msg = "\nmerged_clusters finished succesfully for group " + group_name + "!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nmerged_clusters did not finished succesfully for group "+ group_name + ".\n"
                    log.write(msg)
                    return False
            try:
                outdir_merged_clusters = os.path.join(outdir, "merged_cluster")
                if not os.path.exists(outdir_merged_clusters):
                    os.makedirs(outdir_merged_clusters, exist_ok=True)
                # If you dont want all samples together (maybe you want to group by condition)
                # Please provide a list of lists with the groups you want to make
                merge_cluster_per_group_out = []
                if list(groups.values())[0] == "all":
                    group = list(self.merge_samples.keys())
                    merge_cluster_per_group_out.append(merge_cluster_per_group(group, "merged_cluster"))
                else:
                    # if len(groups) == len(group_names):
                    for group_name, group in groups.items():
                        # group = groups[i]
                        # group_name = group_names[i]
                        log.write("Running merge_clusters for samples " + str(group) + ", members of group " + group_name)
                        if all([sample_id in self.merge_samples.keys() for sample_id in group]):
                            log.write("All samples of group " + group_name + " have been found registered.")
                            merge_cluster_per_group_out.append(merge_cluster_per_group(group, group_name))
                        else:
                            msg = "Not all the sample ids were found in self.merge_samples. Are you sure you are passing sample ids and not sample names?"
                            log.write(msg)
                    else:
                        msg = "The list of groups and the list of the name of the groups should be of equal length."
                        log.write(msg)
                if all(merge_cluster_per_group_out):
                    delattr(self, "merge_samples")
                    msg = "The merged clusters can be now found in self.merge_samples_groups[group_name] (group_name = merged_cluster for all samples together)"
                    log.write(msg)
                    return True
                else:
                    return False
            except KeyboardInterrupt:
                    msg = Bcolors.HEADER + "User interrupted. Finishing merging clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                    print(msg)
                    log.write(msg)

    def set_merge_clusters(self, outdir_merged_clusters, groups):
        with open(self.logfile, "a") as log:
            def set_merge_cluster_per_group(group, group_name):
                # Cluster objects
                merge_samples_lists_clusters = [sample.clusters for sample_id, sample in self.merge_samples.items() if sample_id in group]
                # Cluster number
                merge_samples_clusters = [cluster.cluster_name.split("merged.clusters_")[1] for merge_samples_list_clusters in merge_samples_lists_clusters for cluster in merge_samples_list_clusters]
                
                # Making a dictionary of type cluster_num : [tsv files (one per sample having this cluster)]
                # To create the cluster objects later
                tsvs = { key : set() for key in merge_samples_clusters }
                for sample_id, samples_clusters in self.merge_samples.items():
                    if sample_id in group:
                        for cluster in samples_clusters.clusters:
                            cluster_num = cluster.cluster_name.split("merged.clusters_")[1]
                            tsvs[cluster_num].add(cluster.tsv)

                # The actual dictionary containing a cluster object made out of the different samples' tsv
                merge_clusters = dict.fromkeys(tsvs.keys())
                for cluster_num, tsv in tsvs.items():
                    merge_clusters[cluster_num] = Cluster(cluster_name = (group_name + "_" + str(cluster_num)), tsv = list(tsvs[cluster_num]), logfile = self.logfile)

                self.merge_samples_groups[group_name] = merge_clusters
                self.outdir_merged_clusters_groups[group_name] = outdir_merged_clusters
                log.write("self.merge_samples_groups[" + group_name + "] is now " + str(merge_clusters) + "\n\n\n")
                
            try:
                if not os.path.exists(outdir_merged_clusters):
                    log.write("Directory not found: " + outdir_merged_clusters + " does not exist.")
                    return 2

                set_merge_cluster_per_group_out = []
                log.write("Groups: " + str(groups.values()) + "\n")
                log.write("Group names: " + str(groups.keys()) + "\n")
                if list(groups.keys())[0] == "all":
                    group = list(self.merge_samples.keys())
                    set_merge_cluster_per_group(group, "merged_cluster")
                else:
                    # if len(groups) == len(group_names):
                    for group_name, group in groups.items():
                        # group = groups[i]
                        # group_name = group_names[i]
                        if all([sample_id in self.merge_samples.keys() for sample_id in group]):
                            set_merge_cluster_per_group(group, group_name)
                        else:
                            msg = "Not all the sample ids were found in self.merge_samples. Are you sure you are passing sample ids and not sample names?"
                            log.write(msg)
                    else:
                        msg = "The list of groups and the list of the name of the groups should be of equal length."
                        log.write(msg)

                delattr(self, "merge_samples")

            except KeyboardInterrupt:
                    msg = Bcolors.HEADER + "User interrupted. Finishing merging clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                    print(msg)
                    log.write(msg)

    def map_clusters(self, mode, outdir, gene_gtf, star_index, RAM, out_tmp_dir=None, unique=False, jobs=1):
        print("Running map_clusters with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                self.map_cluster_results = []
                msg = "Mapping clusters.\n"
                log.write(msg)

                if unique:
                    subdirectory = "unique"
                else:
                    subdirectory = "multiple"
                if mode == "merged":
                    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                        for condition, clusters in self.merge_samples_groups.items():
                            for cluster in clusters.values():
                                fastq_dir = os.path.join(outdir, "merged_cluster/")
                                map_outdir = os.path.join(outdir, "map_cluster/", subdirectory)
                                self.map_cluster_results.append(executor.submit(cluster.map_cluster, "Merged", fastq_dir, map_outdir, gene_gtf, star_index, RAM, out_tmp_dir, unique, self.slurm, self.modules))
                else:
                    if mode == "per_sample":
                        samples_dict = self.samples

                        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                            for sample_id, sample in samples_dict.items():
                                for cluster in sample.clusters:
                                    fastq_dir = os.path.join(outdir, "concatenate_lanes/", sample_id)
                                    outdir_sample = os.path.join(outdir, "map_cluster/", subdirectory, sample_id)
                                    self.map_cluster_results.append(executor.submit(cluster.map_cluster, sample_id, fastq_dir, outdir_sample, gene_gtf, star_index, RAM, out_tmp_dir, unique, self.slurm, self.modules))
                    else:
                        msg = "Please specify a mode (merged/per_sample).\n"
                        print(msg)
                        log.write(msg)
                        return 2

                map_cluster_exit_codes = [i.result()[1] for i in self.map_cluster_results]
                map_cluster_allSuccess = all(exit_code == 0 for exit_code in map_cluster_exit_codes)

                if map_cluster_allSuccess:
                    msg = "\nmap_cluster finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nmap_cluster did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing mapping for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)
        
    def TE_counts_clusters(self, mode, outdir, gene_gtf, te_gtf, unique=False, jobs=1):
        print("Running TE_counts with " + str(jobs) + " jobs.\n")
        with open(self.logfile, "a") as log:
            try:
                self.TE_counts_results = []
                msg = "Quantifying TEs.\n"
                log.write(msg)
                
                if unique:
                    subdirectory = "unique"
                else:
                    subdirectory = "multiple"

                log.write(str(self.merge_samples_groups))
                if mode == "merged":
                    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                        for condition, clusters in self.merge_samples_groups.items():
                            for cluster in clusters.values():
                                bam = os.path.join(outdir, "map_cluster/", subdirectory, (cluster.cluster_name + "_Aligned.sortedByCoord.out.bam"))
                                outdir_sample = os.path.join(outdir, "TE_counts/", subdirectory)
                                self.TE_counts_results.append(executor.submit(cluster.TE_count, self.name, "Merged", bam, outdir_sample, gene_gtf, te_gtf, unique, self.slurm, self.modules))
                else:
                    if mode == "per_sample":
                        samples_dict = self.samples

                        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                            for sample_id, sample in samples_dict.items():
                                for cluster in sample.clusters:
                                    bam = os.path.join(outdir, "map_cluster/", subdirectory, sample_id, (cluster.cluster_name + "_Aligned.sortedByCoord.out.bam"))
                                    outdir_sample = os.path.join(outdir, "TE_counts/", subdirectory, sample_id)
                                    self.TE_counts_results.append(executor.submit(cluster.TE_count, self.name, sample_id, bam, outdir_sample, gene_gtf, te_gtf, unique, self.slurm, self.modules))
                            
                    else:
                        msg = "Please specify a mode (merged/per_sample).\n"
                        print(msg)
                        log.write(msg)
                        return 2
               
                TE_counts_exit_codes = [i.result()[1] for i in self.TE_counts_results]
                TE_counts_all_success = all(exit_code == 0 for exit_code in TE_counts_exit_codes)

                if TE_counts_all_success:
                    msg = "\nTE_counts finished succesfully for all samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nTE_counts did not finished succesfully for all samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing TE_counts for all clusters of all samples before closing." + Bcolors.ENDC + "\n" + "\n"
                print(msg)
                log.write(msg)

    def normalize_TE_counts(self, mode, outdir, groups, unique=False, jobs=1):
        print("Running normalize_TE_counts with " + str(jobs) + " jobs.\n")
        def normalize_TE_counts_per_merged_group(group_name, indir, outdir_norm):
            rdata = os.path.join(self.merge_samples_outdir, (self.name + ".rds"))

            if not os.path.exists("merged_samples_norm_scripts"):
                os.makedirs("merge_samples_norm_scripts", exist_ok=True)
            if not os.path.exists(outdir_norm):
                os.makedirs(outdir_norm, exist_ok=True)
            
            cwd = os.path.dirname(os.path.realpath(__file__))
            cmd = ["Rscript", os.path.join(cwd, "r_scripts/normalize_TEexpression.R"), "-m", "merged", "-g", group_name, "-s", ','.join(groups[group_name]), "-o", outdir_norm, "-i", indir, "-r", rdata, "-n", self.name]

            if self.slurm_path != None:
                cmd = ' '.join(cmd)
                
                job_file =  "merge_samplesNorm_scripts/" + self.name + group_name + "_merge_samples_norm.sh"
                job_id = run_job("normalize_TE_counts", job_file, cmd, self.slurm, self.modules)
                
                msg = sucess_submit("normalize_TE_counts", self.name, job_id)
                print(msg)
                log.write(msg)
                
                exit_code = wait_for_job(job_id)

                msg = check_exit_codes("normalize_TE_counts", ("Experiment " + self.name), job_id, exit_code)
                print(msg)
                log.write(msg)

                if exit_code == 0:
                    exit_code = True
            else:
                subprocess.call(cmd)
            self.merge_normalized_outdir = outdir_norm

            if exit_code:
                msg = "\nTE normalization finished succesfully!\n"
                log.write(msg)
                return exit_code
            else:
                msg = "\nTE normalization did not finished succesfully.\n"
                log.write(msg)
                return exit_code
        with open(self.logfile, "a") as log:
            msg = "Normalizing TE counts.\n"
            log.write(msg)
            try:
                if unique:
                    msg = "\nSorry, normalization only available for multiple mapping.\n"
                    log.write(msg)
                    return False
                else:
                    subdirectory = "multiple"

                indir = os.path.join(outdir, "TE_counts", subdirectory)
                outdir_norm = os.path.join(outdir, "TE_counts_normalized", subdirectory)
                
                self.normalized_results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
                    if mode == "merged":
                        for group_name, group in groups.items():
                            self.normalized_results.append(executor.submit(normalize_TE_counts_per_merged_group, group_name, indir, outdir_norm))
                        normalized_exit_codes = [i.result() for i in self.normalized_results]
                        normalized_all_success = all(normalized_exit_codes)
                    else:
                        for sample in list(self.samples.values()):
                            sample_indir = os.path.join(indir, sample.sample_id)
                            sample_outdir = os.path.join(outdir_norm, sample.sample_id)
                            self.normalized_results.append(executor.submit(sample.normalize_TE_counts, sample_indir, sample_outdir))
                            sample.normalized_outdir = sample_outdir
                        normalized_exit_codes = [i.result() for i in self.normalized_results]
                        normalized_all_success = all(exit_code == 0 for exit_code in normalized_exit_codes)

                if normalized_all_success:
                    msg = "\nTE normalization finished succesfully for all groups/samples!\n"
                    log.write(msg)
                    return True
                else:
                    msg = "\nTE normalization did not finished succesfully for all groups/samples.\n"
                    log.write(msg)
                    return False
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                print(msg)
                log.write(msg)

    def process_clusters(self, mode, outdir, gene_gtf, te_gtf, star_index, RAM, groups = {"all" : "merged_cluster"}, out_tmp_dir = None, unique=False, jobs=1, tsv_to_bam = True, filter_UMIs = True, bam_to_fastq = True, concatenate_lanes = True, merge_clusters = True, map_cluster = True, TE_counts = True, normalize_TE_counts = True):
        with open(self.logfile, "a") as log:
            msg = "Running whole pipeline.\n"
            log.write(msg)
            finished_on_the_run_merge_clusters = False
            try:
                if mode == "merged":
                    samples_dict = self.merge_samples
                else:
                    if mode == "per_sample":    
                        samples_dict = self.samples
                    else:
                        msg = "Please specify a mode (merged/per_sample)"
                        print(msg)
                        log.write(msg)
                        return 2

                if tsv_to_bam:
                    current_instruction = "tsv_to_bam"
                    msg = "Running " + current_instruction
                    log.write(msg)
                    tsv_to_bam = self.tsv_to_bam_clusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not tsv_to_bam:
                        msg = "Error in tsv_to_bam"
                        print(msg)
                        log.write(msg)
                        return 1

                if filter_UMIs:
                    current_instruction = "filter_UMIs"
                    msg = "tsv_to_bam finished! Moving on to " + current_instruction
                    log.write(msg)
                    filter_UMIs = self.filter_UMIs_clusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not filter_UMIs:
                        msg = "Error in filter_UMIs"
                        print(msg)
                        log.write(msg)
                        return 1
                    
                if bam_to_fastq:
                    current_instruction = "bam_to_fastq"
                    msg = "filter_UMIs finished! Moving on to " + current_instruction
                    log.write(msg)
                    bam_to_fastq = self.bam_to_fastq_clusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not bam_to_fastq:
                        msg = "Error in bam_to_fastq"
                        print(msg)
                        log.write(msg)
                        return 1

                if concatenate_lanes:
                    current_instruction = "concatenate_lanes"
                    msg = "bam_to_fastq finished! Moving on to " + current_instruction
                    log.write(msg)
                    concatenate_lanes = self.concatenate_lanes_clusters(mode = mode, outdir = outdir, jobs = jobs)
                    if not concatenate_lanes:
                        msg = "Error in concatenate_lanes"
                        print(msg)
                        log.write(msg)
                        return 1

                if mode == "merged" and merge_clusters:
                    self.merge_samples_groups = dict.fromkeys(list(groups.keys()))
                    self.outdir_merged_clusters_groups = dict.fromkeys(list(groups.keys()))

                    current_instruction = "merge_clusters"
                    msg = "concatenate_lanes finished! You selected merged as your mode, so moving on to " + current_instruction
                    log.write(msg)
                    merge_clusters = self.merge_clusters(outdir = outdir, groups = groups)
                    if not merge_clusters:
                        msg = "Error in merge_clusters"
                        print(msg)
                        log.write(msg)
                        return 1
                    finished_on_the_run_merge_clusters = True

                if map_cluster:
                    if mode == "merged" and not finished_on_the_run_merge_clusters:
                        self.merge_samples_groups = dict.fromkeys(list(groups.keys()))
                        self.outdir_merged_clusters_groups = dict.fromkeys(list(groups.keys()))

                        self.set_merge_clusters(os.path.join(outdir, "merged_cluster"), groups = groups)
                        finished_on_the_run_merge_clusters = True

                    current_instruction = "map_cluster"
                    msg = "merge_clusters finished! Moving on to " + current_instruction
                    log.write(msg)
                    map_cluster = self.map_clusters(mode = mode, outdir = outdir, gene_gtf = gene_gtf, star_index = star_index, RAM = RAM, out_tmp_dir = out_tmp_dir, unique = unique, jobs = jobs)
                    if not map_cluster:
                        msg = "Error in map_cluster"
                        print(msg)
                        log.write(msg)
                        return 1

                if TE_counts:
                    if mode == "merged" and not finished_on_the_run_merge_clusters:
                        self.merge_samples_groups = dict.fromkeys(list(groups.keys()))
                        self.outdir_merged_clusters_groups = dict.fromkeys(list(groups.keys()))

                        self.set_merge_clusters(os.path.join(outdir, "merged_cluster"), groups = groups)
                        finished_on_the_run_merge_clusters = True

                    current_instruction = "TE_counts"
                    msg = "map_cluster finished! Moving on to " + current_instruction
                    log.write(msg)
                    TE_counts = self.TE_counts_clusters(mode = mode, outdir = outdir, gene_gtf = gene_gtf, te_gtf = te_gtf, unique = unique, jobs = jobs)
                    if not TE_counts:
                        msg = "Error in TE_counts"
                        print(msg)
                        log.write(msg)
                        return 1

                if normalize_TE_counts:
                    if mode == "merged" and not finished_on_the_run_merge_clusters:
                        self.merge_samples_groups = dict.fromkeys(list(groups.keys()))
                        self.outdir_merged_clusters_groups = dict.fromkeys(list(groups.keys()))

                        self.set_merge_clusters(os.path.join(outdir, "merged_cluster"), groups = groups)
                        finished_on_the_run_merge_clusters = True

                    current_instruction = "normalize_TE_counts"
                    msg = "TE_counts finished! Moving on to " + current_instruction
                    log.write(msg)
                    normalize_TE_counts = self.normalize_TE_counts(mode = mode, outdir = outdir, groups = groups, unique = unique, jobs = jobs)
                    if not normalize_TE_counts:
                        msg = "Error in normalize_TE_counts"
                        print(msg)
                        log.write(msg)
                        return 1

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted. Finishing instruction " + current_instruction + " for all clusters of all samples before closing." + Bcolors.ENDC + "\n"
                print(msg)
                log.write(msg)

            
