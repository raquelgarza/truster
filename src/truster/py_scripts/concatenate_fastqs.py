#!/usr/bin/python

import sys
import os
import gzip
import getopt
from Bio import SeqIO
from itertools import count

def main(argv):
    fastqs_in = ''
    fastq_out = ''
    sample_id = ''
    cluster_name = ''
    library_names = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:c:l:",["fastqs_in=","fastq_out=","sample_id=","cluster_name=","library_names="])
    except getopt.GetoptError:
        print("%s: %s\n" % (args[0], e.msg))
        print('concatenate_fastqs.py -i <fastqs_in> -o <fastq_out> -s <sample_id> -c <cluster_name> -l <library_names>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('concatenate_fastqs.py -i <fastqs_in> -o <fastq_out> -s <sample_id> -c <cluster_name> -l <library_names>')
            sys.exit()
        elif opt in ("-i", "--fastqs_in"):
            fastqs_in = arg
        elif opt in ("-o", "--fastq_out"):
            fastq_out = arg
        elif opt in ("-s", "--sample_id"):
            sample_id = arg
        elif opt in ("-c", "--cluster_name"):
            cluster_name = arg
        elif opt in ("-l", "--library_names"):
            library_names = arg

    fastqs_in = fastqs_in.split(",")
    library_names = library_names.split(",")
    
    if not os.path.exists(os.path.dirname(fastq_out)):
    	os.makedirs(os.path.dirname(fastq_out))

    with gzip.open(fastq_out, "at") as fq_out:
        for i in range(0,len(fastqs_in)):
            fastq = fastqs_in[i]
            library_name = library_names[i]
            c = count()
            with gzip.open(fastq, "rt") as fq_in:
                for record in SeqIO.parse(fq_in, "fastq"):
                    record.id = f"{sample_id}_{cluster_name}_{library_name}_readnum{next(c)}"
                    SeqIO.write(record, fq_out, "fastq")

if __name__ == "__main__":
   main(sys.argv[1:])
