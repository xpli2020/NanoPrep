import os
import re
import pandas as pd
import subprocess
from collections import defaultdict
import shutil
from Bio import SeqIO

def searchSubsLevels(dir):
    dir =re.sub("/$", "", dir)

    each_level = []
    for i, _, _ in os.walk(dir):
        each_level.append(i)
    
    level_folders = defaultdict(list)
    level_folders_path = defaultdict(list)
    for i in each_level:

        level_folders[len(i.split("/"))-1].append(i.split("/")[len(i.split("/"))-1])
        level_folders_path[len(i.split("/"))-1].append(i)
        
    return level_folders, level_folders_path


def getFiles(out_searchSublevels, levels = 1):
    """[summary]

    Args:
        out_searchSublevels ([type]): function searchSubsLevels output

    Returns:
        [type]: a dictionary {sublevel folder: files}
    """

    folders = out_searchSublevels[0][levels]
    paths = out_searchSublevels[1][levels]

    folder_file_dict = defaultdict(list)
    for i,v in enumerate(folders):
        searchfile = os.listdir(paths[i])
        for j in searchfile:
            filepath = os.path.join(paths[i], j)
            if os.path.isfile(filepath):
                folder_file_dict[v].append(filepath)

    return folder_file_dict


def recursive_getfiles(parentpath, subdir = True):
    # max level = 2
    
    if not os.path.exists(parentpath):
        return "Path not found"
    
    p_path = os.path.abspath(parentpath)

    all_files = {}
    subs = os.listdir(p_path)
    if subdir:
        for i in subs:
            i_path = os.path.join(p_path, i)
            d1, subfolder, _ = next(os.walk(i_path))
            subfolder_path = os.path.join(d1, subfolder[0])
            _,_,q = next(os.walk(subfolder_path))
            # all_files[i] = [ os.path.join(subfolder_path, k) for k in q if q != "none.fastq"]
            all_files[i] = [ os.path.join(subfolder_path, k) for k in q ]
            
    return all_files



"""Preprocess the nanopore reads

need a input file mapping using fq_src_mapping.py

1. combined reads
2. demultiplex and trimming with porechop or qcat
3. nanoplot
4. nanoFilt
"""

def nano_map2dict(mapping, sep="\t", header = None):
    """turn mapping file into dictionary

    Args:
        mapping (txt, csv): [description]
    """

    df = pd.read_csv(mapping, sep=sep, names=["Sample", "Read_path"], header=header)

    read_mapping = defaultdict(list)
    for _, content in df.iterrows():
        read_mapping[content["Sample"]].append(content["Read_path"])

    return read_mapping

def nano_groupreads(mapping, OUTPUT="Grouped_reads", groupfile=None, keepSubdir=False, changeName=True, changePrefix="", exist_OK = False):
    """This function is used to group reads, i.e 16S and ITS by the groupfile, or change folder and read names if necessary

    Args:
        mapping ([type]): [description]
        OUTPUT (str, optional): [description]. Defaults to "grouped_reads".
        groupfile ([type], optional): [provide a tab file to specify groups, need "GROUP" column and "SampleID" should match to the mapping input keys]. Defaults to None
        keepSubdir (Boolean): if keep the orginial subfolders, or put sub files together by group
        changeName: change folder name and file name
        changePrefix: add prefix to folder name and file name
    """
    # move files first -> then classify to groups

    os.makedirs(OUTPUT, exist_ok=exist_OK)
    root_dir = os.path.abspath(OUTPUT)


    # add summary report
    raw_read = defaultdict(list)

    moved_reads = nano_copysource(mapping, OUTPUT=OUTPUT, exist_ok=True)

    # moved_reads : {"path/barcode01":[...]}

    groupGuide = pd.read_csv(groupfile, sep="\t", names=["SampleID", "Barcodes", "Group"])
    # get unique groups
    # {16S: [barcode01, barcode02, ...]}
    GRPs = defaultdict(list)
    for row, content in groupGuide.iterrows():
        GRPs[content["Group"]].append(str(content["Barcodes"])+"_"+str(content["SampleID"]))

    
    # make subfolders by group
    for k, v in GRPs.items():
        # k: group i.e. 16S or ITS
        # v: barcode_sampleID
        if not os.path.exists(k):
            # make group subfoler
            sam = os.path.join(root_dir, os.path.basename(k))
            os.makedirs(sam)
            
        # extract barcode and sampleID
        for j in v:
            
            b_s = j.split("_")
            bar = b_s[0]
            sampleid = b_s[1]
            
            # match with moved_reads
            path_moved = os.path.join(root_dir, bar)
            reads_need_moving = moved_reads[path_moved] # a list of files with absolute path
            
            for mread in reads_need_moving:
            
                # if keepsubdir
                print("\n......keeping sub folders ......\n")
                if keepSubdir:
                    # make sub folder under each group folder
                    if changeName: # change subfolder name
                        print("\n...rename sub folders and files...\n")
                        new_name =changePrefix+sampleid+"_"+bar
                        sub_sam = os.path.join(sam, new_name)
                        chgn_read = os.path.join(sub_sam, new_name+"."+os.path.basename(mread).split(".")[1])
                    else:
                        sub_sam = os.path.join(sam, bar)
                        chgn_read = os.path.join(sub_sam, os.path.basename(mread))

                    if not os.path.exists(sub_sam): # create sub folders
                        os.makedirs(sub_sam, exist_ok=exist_OK)
                        
                        
                    # now move the reads
                    raw_read[sub_sam].append(chgn_read)
                    print("\n.Grouping.\n")
                    command = " ".join(["cp", mread, chgn_read])
                    subprocess.run(command, shell=True)
                    
                    

                else: # not keepsubdir
                    print("\n......Not keeping sub folders ......\n")
                    if changeName: # change subfolder name
                        new_name =changePrefix+sampleid+"_"+bar
                        chgn_read = os.path.join(sam, new_name+"."+os.path.basename(mread).split(".")[1])
                    else:
                        chgn_read = os.path.join(sam, os.path.basename(mread))
                                        # now move the reads
                            
                    raw_read[sam].append(chgn_read)
                    print("\n.Grouping.\n")
                    command = " ".join(["cp", mread, chgn_read])
                    subprocess.run(command, shell=True)


    # remove the copied file to save space
    print("!!remove copied reads to save space\n")
    for i in moved_reads.keys():
        command = " ".join(["rm", i, "-rf"])
        subprocess.run(command, shell=True)


    return raw_read

def nano_copysource(mapping, OUTPUT="copied_fastq_reads", exist_ok=False):
    """Copy source fastq files to current working direction, to avoid mess up original files

    mapping: mapping file (Sample, path fastq)
    return a dcitionary contains {uplevel folder name path: new path fastq}

    """

    os.makedirs(OUTPUT, exist_ok=exist_ok)
    root_dir = os.path.abspath(OUTPUT)

    # add summary report
    raw_read = defaultdict(list)

    if type(mapping) == defaultdict:

        fastq_source = mapping

        for k, v in fastq_source.items():
            for j in v:
                sam = os.path.join(root_dir, os.path.basename(k))
                read_from = j

                if os.path.exists(read_from):
                    read_to = os.path.join(sam, os.path.basename(read_from))
                    print(f"\nCopying  sample {read_from} to working directory {sam}")
                else:
                    print("read path does not exist")
                    return "read path does not exist"

                if not os.path.exists(sam):
                    os.makedirs(sam)

                raw_read[sam].append(read_to)
                command = " ".join(["cp", read_from, read_to])
                subprocess.run(command, shell=True)
            

    else:
        fastq_source = pd.read_csv(mapping, sep="\t", names=["Sample", "Read_path"], header=None)
        print(f"Total reads from source {fastq_source.shape[0]}")

        for row, content in fastq_source.iterrows():

            sam = os.path.join(root_dir, os.path.basename(content["Sample"]))
            read_from = content["Read_path"]
            if os.path.exists(read_from):
                read_to = os.path.join(sam, os.path.basename(read_from))
                print(f"\n{row} Copying  sample {read_from} to working directory {sam}")
            else:
                print("read path does not exist")
                return "read path does not exist"

            if not os.path.exists(sam):
                os.makedirs(sam)

            raw_read[sam].append(read_to)
            command = " ".join(["cp", read_from, read_to])
            subprocess.run(command, shell=True)
            
    print("Completed")
    return raw_read


def nano_concat(sample_read_dict, mergeby="Sample", output="Concat_reads", filetype="fastq", exist_ok=False):
    """[summary] concat reads, if by sample, concat reads under each sample, if by all, concat all reads 

    Args:
        sample_read_dict (dict): [sample read dictionary] 
        mergeby (str, optional): [description]. Defaults to "Sample".
    """

    assert type(sample_read_dict) == defaultdict, "input is a python dictionary with {sample:path, ...}"

    try:
        os.makedirs(output, exist_ok=exist_ok)
    except:
        print(f"{output} exist")
        return


    out = os.path.abspath(output)

    read_mapping = defaultdict(list)
    if mergeby == "Sample":
        print("Concat reads by sample...")
        for sample, reads in sample_read_dict.items():

            sample_path = os.path.join(out, os.path.basename(sample))
            os.makedirs(sample_path, exist_ok=True)
            outread = "concat_"+os.path.basename(sample)+"."+filetype
            outread_path = os.path.join(sample_path, outread)
            read_mapping[sample_path].append(outread_path)
            # resolve arguments too long eror with cat
            uplevelfrom = os.path.dirname(reads[0]) 
            print(f"//Find files from source: {uplevelfrom}")
            command = " ".join(["cat", f"{uplevelfrom}/*.{filetype}", ">", outread_path])
            print(command)
            print("\n")
            subprocess.run(command, shell=True)

    elif mergeby == "All":
        print("Concat all reads...")
        all_reads = [ j for i in sample_read_dict.values() for j in i]
        all_reads_upper = { os.path.dirname(i)+f"/*.{filetype}" for i in all_reads }
        outread = "all_concat_reads"+"."+filetype
        outread_path = os.path.join(out, outread)
        read_mapping[out].append(outread_path)
        # resolve arguments too long error with cat
        command = " ".join(["cat", *all_reads_upper, ">", outread_path])
        print(command)
        print("\n")
        subprocess.run(command, shell=True)

    print("\nFinished concat!")
    return read_mapping


def nano_runProg(sample_read_dict, *args, method="qcat", inputflag="-f", outputflag="-o",outfolder=None, outfolder_suffix="out", forceout = False, stdin = False, xlabel="", suffix_change=False, outfmt="", exist_ok=True, switchinput_output=False, catchstdout=False, stdout_fmt=""):
    """Demultiplex with qcat or porechop

    Args:
        sample_read_dict ([type]): [description]
        method (str, optional): [description]. Defaults to "qcat".
        forceout: if the program do not create out file, you can set it to True
        stdin: standard in format: load options first, default false
        suffix_change: change file type
        outfmt: output format to change to
        exist_ok: if override existing folder
        switchinput_output: some program where stdin is the last command
        xlabel: change name on the output file
        catchstdout: if to catch stdout in a file, stdout file name will be the same as output name, except for the suffix, you can change it with stdout_fmt
        stdout_fmt: what file type you want to store stdout

    Returns:
        defaultdict list: with folder:[files]
    """

    
#     if method not in ["qcat", "porechop"]:
#         print(f"{method} not supported currently")
#         return
    
    # check if install programs

    results = defaultdict(list)

    if shutil.which(method) is None:
        print(f"Cannot find {method} in the system path. Please install the program")
        return

    else:

        if outfolder is None:
            out = f"{method}_{outfolder_suffix}"
        else:
            out = f"{outfolder}_{outfolder_suffix}"

        os.makedirs(out, exist_ok=exist_ok)
        out_path = os.path.abspath(out)
        
        for sample, reads in sample_read_dict.items():
            
            sample_out = os.path.join(out_path, os.path.basename(sample))
            
            if forceout:
                os.makedirs(sample_out, exist_ok=True)
            
            for read in reads:

                if stdin:

                    basecommand = [method]
                    basecommand.append(*args)
                    # deal with suffix, file type
                    if suffix_change:
                        read_label = xlabel + "".join(os.path.basename(read).split(".")[:-1]) +outfmt
                    else:
                        read_label = xlabel + os.path.basename(read)

                    read_label_path = os.path.join(sample_out, read_label)

                    if switchinput_output:
                        basecommand.extend([outputflag, read_label_path, inputflag, read])
                    else:
                        basecommand.extend([inputflag, read, outputflag, read_label_path])
                    
                    new_command = " ".join(basecommand)
                    
                    print(f"\nRunning program: {method} on: {sample}")
                    print(f">>Input: {read}")
                    print(f">>Output: {sample_out}")
                    print(f"[ Command: {new_command} ]")

                    if catchstdout:
                        
                        stdout_name = xlabel + "".join(os.path.basename(read).split(".")[:-1]) + "." + stdout_fmt
                        stdout_file_path = os.path.join(sample_out, stdout_name)
                        print(f"Catching stdout in {stdout_file_path}")
                        with open(stdout_file_path, "w") as stdoutf:
                            Result = subprocess.run(new_command, shell=True, stdout=stdoutf, stderr=subprocess.PIPE)

                    else:
                        Result = subprocess.run(new_command, shell=True, stderr=subprocess.PIPE)

                    print(Result.stderr)

                else:
                    # when the program generates folder for you: sample out
                    basecommand = [method, inputflag, read, outputflag, sample_out]
                    basecommand.append(*args)
                    
                    new_command = " ".join(basecommand)
                    
                    print(f"\nRunning program: {method} on: {sample}")
                    print(f">>Input: {read}")
                    print(f">>Output: {sample_out}")
                    print(f"[ Command: {new_command} ]")
            
                    Result = subprocess.run(new_command, shell=True, stderr=subprocess.PIPE)
                    print(Result.stderr)

            outlist = os.listdir(sample_out)
            for i in outlist:
                results[sample_out].append(os.path.join(sample_out, i))


    return results

# filter barcode according to the guidefile

def regexBarcode(barcode_file_name, regpattern="barcode([0-9]+)", groupN = 1):
    """Uniform output barcode name to BC1, BC2

    Args:
        barcode_file_name ([type]): [description]
        method (str, optional): [description]. Defaults to "qcat".

    Returns:
        [type]: [description]
    """
    # BC01, BC11
    # BC01 -> BC1

    regex = re.compile(regpattern)

    if regex.search(barcode_file_name) is not None:
        label = int(regex.search(barcode_file_name).group(groupN))
        newlabel = "BC"+str(label)
        return newlabel
    else:
        return ""



def nano_filterBC(sample_read_dict, guide_file, copy = True, sampleCOL="Samples_folder", barcodeCOL = "Barcodes", delim="\t", copyout = "Filtered_demutiplex_barcode", mappingout = "filtered_barcode.mapping", exist_ok=False):
    """Filter barcode by guided file, returns a folder contains only barcodes in the guide file

    Args:
        sample_read_dict ([type]): [description]
        guide_file ([type]): [description]
        copy (bool): if copy filtered barcode file into a new folder
        sampelCOL (str): the column in the guide file pointing to the sample folder
        barcodeCOL (str): the column in the guide file pointing to the barcode
        delim (str): the delimiter in the guid file
    """

    guide_df = pd.read_csv(guide_file, sep=delim)

    guided_sam_bc = defaultdict(list)

    for _, content in guide_df.iterrows():

        barcode = content[barcodeCOL]
        sample = content[sampleCOL]
        guided_sam_bc[sample].append(barcode)



    filtered_sample_read_dict = defaultdict(list)

    for key, value in sample_read_dict.items():

        sample2 = os.path.basename(key)
        
        for bar in value:
            # convert barcode name to match guided file, none.fastq -> ""
            barlabel = regexBarcode(bar)
            if barlabel in guided_sam_bc[sample2]:
                filtered_sample_read_dict[sample2].append(bar)


    if copy:
        
        os.makedirs(copyout, exist_ok=exist_ok)
        copyoutpath = os.path.abspath(copyout)
        print(f"Copy filterd files into folder {copyoutpath} ...")
        new_filtered_dict = defaultdict(list)

        for key, value in filtered_sample_read_dict.items():
            # make sub directory
            subsam_path = os.path.join(copyoutpath, key)
            os.makedirs(subsam_path)
            print(f"Working on {key}")

            for bar in value:
                barlabel = regexBarcode(bar)
                barout = os.path.join(subsam_path, os.path.basename(bar))
                print(f">> Copying barcode {barlabel}: from {bar} to {barout}")
                command = " ".join(["cp", bar, barout])
                subprocess.run(command, shell=True)

                new_filtered_dict[subsam_path].append(barout)

        with open(mappingout, "w") as f:
            print("Samples_folder", "Barcodes", sep="\t", file=f)
            for k, v in new_filtered_dict.items():
                for j in v:
                    print(k, j, sep="\t", file=f)

        return new_filtered_dict

    else:

        with open(mappingout, "w") as f:
            print("Samples_folder", "Barcodes", sep="\t", file=f)
            for k, v in filtered_sample_read_dict.items():
                for j in v:
                    print(k, j, sep="\t", file=f)

        return filtered_sample_read_dict


        
#################################################################

# nano_fish() - separate mixed libraries?

def process(lines=None):
    ks = ["name", "sequence", "optional", "quality"]
    return {k:v for k, v in zip(ks, lines)}

def changeFastqHeader(fastq, newrecordname=None, n=4, out="changeFastqHeader.fastq"):

    with open(out, "w") as out:
        with open(fastq, "r") as fh:
            lines = []
            print(f"Opening {fastq}")
            for line in fh:
                lines.append(line.rstrip())
                if len(lines) == n:
                    record = process(lines)
                    if newrecordname is not None:
                        record_id = record["name"].split()[0]
                        newheader = newrecordname + "|" + record_id
                        record["name"] = newheader

                    print("@"+record["name"], sep="", file=out)
                    print(record["sequence"], sep="", file=out)
                    print(record["optional"], sep="", file=out)
                    print(record["quality"], sep="", file=out)
                    lines = []



def nano_modiheader(sample_read_dict, filetype="fasta", out="Modified_header", exist_ok=True):
    """Modify read name to SampleID|BC|readid

    Args:
        sample_read_dict ([type]): [description]
        filetype (str, optional): [description]. Defaults to "fasta".
        out (str, optional): [description]. Defaults to "Modified_header".
        exist_ok (bool, optional): [description]. Defaults to True.
    """


    assert type(sample_read_dict) == defaultdict, "input is a python dictionary with {sample:path, ...}"

    os.makedirs(out, exist_ok=exist_ok)

    outpath = os.path.abspath(out)

    read_mapping = defaultdict(list)


    # reads and sample
    # load each fasta file
    # each read would have an ID >sample_folder|barcode|readID
    # here barcode would be converted to barcode 1 to BC1

    for sample, reads in sample_read_dict.items():
        sample_folder = os.path.basename(sample)

        sample_folder_path = os.path.join(outpath, sample_folder)

        if os.path.exists(sample_folder_path):
            print("Sample exist, will overide")

        os.makedirs(sample_folder_path, exist_ok=True)

        for read in reads:
            barcode = os.path.basename(read)
            
            newbarcode = regexBarcode(barcode)
            header1 = "|".join([sample_folder, newbarcode])
            outputname = sample_folder + "_" + newbarcode + f".{filetype}"
            barcode_path = os.path.join(sample_folder_path, outputname)

            read_mapping[sample_folder_path].append(barcode_path)

            if filetype in ["fasta","fa", "fna"]:
                seqRecord = SeqIO.parse(read, format= filetype) # this is a generator
                with open(barcode_path, "w") as outf:
                    for record in seqRecord:
                        record.id = header1 + "|" + record.id
                        seqs = record.seq 
                        print(">"+record.id, sep="", file=outf)
                        print(seqs, sep="", file=outf)

            elif filetype in ["fastq", "fq"]:
                changeFastqHeader(read, newrecordname=header1, out=barcode_path)

    return read_mapping



def nano_krakenReport(sample_read_dict, saveSummary="classfication_summary.txt"):
    """A summary of classfied and unclassifed reads all together
    """
    
    classified = "C"
    unclassified = "U"

    classified_count = 0
    unclassified_count = 0
    total_count = 0


    with open(saveSummary, "w") as su:
        print("Sample", "Kraken_file", "Total", "Classified", "Unclassified", "Classified_perc", "Unclassified_perc", sep="\t", file=su)
        for sample, files in sample_read_dict.items():
            for f in files:
                df = pd.read_csv(f, sep="\t", header=None, names=["V0", "V1", "V2", "V3", "V4"])

                n_classified =  sum(df["V0"] == classified)
                n_unclassified = sum(df["V0"] == unclassified)
                total_n = df.shape[0]

                classified_count += n_classified
                unclassified_count += n_unclassified
                total_count += total_n
                print(f"{os.path.basename(sample)}", f"{os.path.basename(f)}", f"{total_n}", f"{n_classified}", f"{n_unclassified}", f"{100*n_classified/total_n:.4f}", f"{100*n_unclassified/total_n:.4f}",sep="\t", file=su)


        print(f"Total", f"all_files", f"{total_count}", f"{classified_count}", f"{unclassified_count}", f"{100*classified_count/total_count:.4f}", f"{100*unclassified_count/total_count:.4f}",sep="\t", file=su)
