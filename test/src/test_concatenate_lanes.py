#!/bin/python
# arrange
import truster
import os
import gzip

raw_path = "/Volumes/My Passport/trusTEr/test/raw/"
lunarc = "config_files/lunarc_config.json"
modules = "config_files/software_modules.json"

test = truster.Experiment("test", lunarc, modules)
test.register_samples_from_path(raw_path)

merge_samples_outdir = "/Volumes/My Passport/trusTEr/test/3_merge_samples/"
test.set_merge_samples_outdir(merge_samples_outdir)

test.slurm = None
test.modules = None
test.concatenate_lanes_clusters(mode="merged", outdir="/Volumes/My Passport/trusTEr/test/3_merge_samples/cluster_pipeline/")

# Expected
expected_merge_clusters_outdir = "/Volumes/My Passport/trusTEr/test/3_merge_samples/cluster_pipeline/concatenate_lanes/expected_results"
with gzip.open(os.path.join(expected_merge_clusters_outdir, "sample_A", "sample_A_merged.clusters_0_R2.fastq.gz"), "rt") as sample_A_0:
    expected_sample_A_0_file = sample_A_0.readlines()
# Actual
merge_clusters_outdir = "/Volumes/My Passport/trusTEr/test/3_merge_samples/cluster_pipeline/concatenate_lanes"
with gzip.open(os.path.join(merge_clusters_outdir, "sample_A", "sample_A_merged.clusters_0_R2.fastq.gz"), "rt") as sample_A_0:
    sample_A_0_file = sample_A_0.readlines()

# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "sample_A", "sample_A_merged.clusters_1_R2.fastq.gz"), "rt") as sample_A_1:
    expected_sample_A_1_file = sample_A_1.readlines() 
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "sample_A", "sample_A_merged.clusters_1_R2.fastq.gz"), "rt") as sample_A_1:
    sample_A_1_file = sample_A_1.readlines()

# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "sample_B", "sample_B_merged.clusters_0_R2.fastq.gz"), "rt") as sample_B_0:
    expected_sample_B_0_file = sample_B_0.readlines()
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "sample_B", "sample_B_merged.clusters_0_R2.fastq.gz"), "rt") as sample_B_0:
    sample_B_0_file = sample_B_0.readlines()

# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "sample_B", "sample_B_merged.clusters_1_R2.fastq.gz"), "rt") as sample_B_1:
    expected_sample_B_1_file = sample_B_1.readlines() 
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "sample_B", "sample_B_merged.clusters_1_R2.fastq.gz"), "rt") as sample_B_1:
    sample_B_1_file = sample_B_1.readlines()

# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "sample_C", "sample_C_merged.clusters_0_R2.fastq.gz"), "rt") as sample_C_0:
    expected_sample_C_0_file = sample_C_0.readlines() 
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "sample_C", "sample_C_merged.clusters_0_R2.fastq.gz"), "rt") as sample_C_0:
    sample_C_0_file = sample_C_0.readlines()

# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "sample_C", "sample_C_merged.clusters_1_R2.fastq.gz"), "rt") as sample_C_1:
    expected_sample_C_1_file = sample_C_1.readlines()
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "sample_C", "sample_C_merged.clusters_1_R2.fastq.gz"), "rt") as sample_C_1:
    sample_C_1_file = sample_C_1.readlines()

# assert
print(sorted(expected_sample_A_0_file) == sorted(sample_A_0_file))
print(sorted(expected_sample_A_1_file) == sorted(sample_A_1_file))
print(sorted(expected_sample_B_0_file) == sorted(sample_B_0_file))
print(sorted(expected_sample_B_1_file) == sorted(sample_B_1_file))
print(sorted(expected_sample_C_0_file) == sorted(sample_C_0_file))
print(sorted(expected_sample_C_1_file) == sorted(sample_C_1_file))




