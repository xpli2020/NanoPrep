#!/usr/bin/env python3

import os
import pandas as pd
import argparse
from collections import defaultdict



def makemapfile(fq_file_directory, outfile="mapping.txt", filter=False, filtertype="fastq"):

    
    fq_file_directory = os.path.abspath(fq_file_directory)
    newfile_dict = defaultdict(list)

    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(fq_file_directory):
        listOfFiles += [os.path.join(dirpath, file) for file in filenames]

    for i in listOfFiles:
        dirname = os.path.dirname(i)
        newfile_dict[dirname].append(i)


    with open(outfile, "w") as outf:
        for k, v_list in newfile_dict.items():
            for v in v_list:
                if filter:
                    if v.split(".")[-1] == filtertype:
                        print(k, v, sep="\t", file=outf)
                else:
                    print(k, v, sep="\t", file=outf)

    return newfile_dict

    
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-fqdir', dest="fqdir", help="fq file directory")
    parser.add_argument('-out', dest="out", default="mapping.txt", help="mapping file name out")
    parser.add_argument('-filter', dest="filter", action="store_true", help="if filter file")
    parser.add_argument('-ftype', dest="ftype", default="fastq", help="file type")

    args = parser.parse_args()

    fqdir = args.fqdir
    out = args.out
    fil = args.filter
    ftype = args.ftype
    
    df = makemapfile(fqdir, out, fil, ftype)



