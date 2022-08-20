#!/bin/python
# arrange
import truster
import os
import gzip

raw_path = "/Volumes/MyPassport/trusTEr/test/raw/"
lunarc = "config_files/lunarc_config.json"
modules = "config_files/software_modules.json"

test = truster.Experiment("test", lunarc, modules)
test.register_samples_from_path(raw_path)

merge_samples_outdir = "/Volumes/MyPassport/trusTEr/test/3_merge_samples/"
test.set_merge_samples_outdir(merge_samples_outdir)

merge_clusters_outdir = "/Volumes/MyPassport/trusTEr/test/3_merge_samples/cluster_pipeline/"
groups = {"group_A" : ["sample_A", "sample_B"], "group_B" : ["sample_C"]}

test.merge_samples_groups = dict.fromkeys(list(groups.keys()))
test.outdir_merged_clusters_groups = dict.fromkeys(list(groups.keys()))
# act
test.merge_clusters(outdir = merge_clusters_outdir, groups = groups) 

# Check group A, cluster 0
# Expected
expected_merge_clusters_outdir = "/Volumes/MyPassport/trusTEr/test/3_merge_samples/cluster_pipeline/expected_results/merged_cluster/"
with gzip.open(os.path.join(expected_merge_clusters_outdir, "group_A_0_R2.fastq.gz"), "rt") as group_A_0:
    expected_group_A_0_file = group_A_0.readlines()
# Actual
merge_clusters_outdir = "/Volumes/MyPassport/trusTEr/test/3_merge_samples/cluster_pipeline/merged_cluster"
with gzip.open(os.path.join(merge_clusters_outdir, "group_A_0_R2.fastq.gz"), "rt") as group_A_0:
    group_A_0_file = group_A_0.readlines()

# Check group A, cluster 1
# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "group_A_1_R2.fastq.gz"), "rt") as group_A_1:
    expected_group_A_1_file = group_A_1.readlines()
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "group_A_1_R2.fastq.gz"), "rt") as group_A_1:
    group_A_1_file = group_A_1.readlines()

# Check group B, cluster 0
# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "group_B_0_R2.fastq.gz"), "rt") as group_B_0:
    expected_group_B_0_file = group_B_0.readlines()
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "group_B_0_R2.fastq.gz"), "rt") as group_B_0:
    group_B_0_file = group_B_0.readlines()

# Check group B, cluster 1
# Expected
with gzip.open(os.path.join(expected_merge_clusters_outdir, "group_B_1_R2.fastq.gz"), "rt") as group_B_1:
    expected_group_B_1_file = group_B_1.readlines()
# Actual
with gzip.open(os.path.join(merge_clusters_outdir, "group_B_1_R2.fastq.gz"), "rt") as group_B_1:
    group_B_1_file = group_B_1.readlines()

# assert
print(sorted(expected_group_A_0_file) == sorted(group_A_0_file))
print(sorted(expected_group_A_1_file) == sorted(group_A_1_file))
print(sorted(expected_group_A_0_file) == sorted(group_A_0_file))
print(sorted(expected_group_B_1_file) == sorted(group_B_1_file))




