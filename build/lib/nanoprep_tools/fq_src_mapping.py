#!/usr/bin/env python3

import os
import pandas as pd
import argparse
from collections import defaultdict



def makemapfile(fq_file_directory, subdir_contains_fq="fastq_pass", out="mapping.txt"):

    subfolders = []
    reads = []
    fq_file_directory = os.path.abspath(fq_file_directory)
    
    sub = False
    for i, _, _ in os.walk(fq_file_directory):
        if subdir_contains_fq in i:
            sub = True
            updir = os.path.basename(os.path.dirname(i))            
            for j in os.listdir(i):
                if j.endswith(("fq", "fastq", "fastq.gz")):
                    subfolders.append(updir)
                    reads.append(os.path.join(i, j))

                if j.endswith(("fa", "fasta", "fasta.gz", "fna", "fna.gz", "fa.gz")):
                    subfolders.append(updir)
                    reads.append(os.path.join(i, j))
    if not sub:
        sub_folders = os.listdir(fq_file_directory)
        
        for i in sub_folders:
        
            if os.path.isdir(i):
                sub_reads = os.listdir(os.path.join(fq_file_directory, i))
                sample_name = i
                for j in sub_reads:

                    if j.endswith(("fq", "fastq", "fastq.gz")):
                        subfolders.append(sample_name)
                        reads.append(os.path.join(fq_file_directory, i, j))

                    if j.endswith(("fa", "fasta", "fasta.gz", "fna", "fna.gz", "fa.gz")):
                        subfolders.append(sample_name)
                        reads.append(os.path.join(fq_file_directory, i, j))
            else:
                subfolders.append(fq_file_directory)
                if i.endswith(("fq", "fastq", "fastq.gz")):
                    reads.append(os.path.join(fq_file_directory, i))
                if i.endswith(("fa", "fasta", "fasta.gz", "fna", "fna.gz", "fa.gz")):
                    reads.append(os.path.join(fq_file_directory, i))                  

    df = pd.DataFrame({"Sample":subfolders, "Read_path":reads})

    df.to_csv(out, sep="\t", header=None, index=None)

    return df

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-fqdir', dest="fqdir", help="fq file directory")
    parser.add_argument('-out', dest="out", default="mapping.txt", help="mapping file name out")
    # parser.add_argument('-suffix', dest="suffix", default="fastq")
    parser.add_argument('-subdir', dest="subdir", default="fastq_pass", help="subfolder name that contains fastq reads")
    parser.add_argument('--showdf', dest="showdf", action="store_true", help="if want to show df in standardou")

    args = parser.parse_args()

    fqdir = args.fqdir
    out = args.out
    # suf = args.suffix
    subdir = args.subdir
    showdf = args.showdf

    if showdf:
        df = makemapfile(fqdir, subdir, out)
        print(df.head())

    else:
        makemapfile(fqdir, subdir, out)


