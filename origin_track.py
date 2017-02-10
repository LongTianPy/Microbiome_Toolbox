#!/usr/bin/python3
"""
This function tracks if the sequences of one bacteria taxon found are consistent between samples.
That is, to check if the bacteria x of Day 1 the same as that of Day 2.
"""

# IMPORT
import argparse
import sys
import os
from Bio import SeqIO
from os.path import isdir,isfile, join
from os import listdir
import os

# FUNCTIONS
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="Track and compare the origin of a bacteria"
    )
    parser.add_argument("-i", dest="input", help="Directory of the QIIME workspace")
    parser.add_argument("-l", dest="rank", help="Taxonomic rank")
    parser.add_argument("-t", dest="taxon", help="A taxon you want to look at")
    parser.add_argument("-m", dest="mapping", help="Mapping file")
    args = parser.parse_args()
    return args

def taxonomy_rank(rank):
    rank_lower = rank.lower()
    taxonomy = {"kingdom":"D_0__","phylum":"D_1__","class":"D_2__",
                "order":"D_3__","family":"D_4__","genus":"D_5__",
                "strain":"D_6__"}
    return taxonomy[rank_lower]

def extract_sequences(rank, taxon):
    if not isdir("origin_track/"):
        os.mkdir("origin_track/")
    f = open("step1_otus/step1_rep_set.fna", "r")
    step1_otus = list(SeqIO.parse(f, "fasta"))
    f.close()
    f = open("step4_otus/step4_rep_set.fna", "r")
    step4_otus = list(SeqIO.parse(f, "fasta"))
    f.close()
    otus = step1_otus + step4_otus
    otu_dict = {i.id: [i.description, str(i.seq)] for i in otus}
    f = open("uclust_assigned_taxonomy/rep_set_tax_assignments.txt", "r")
    taxa_table = [i.strip().split("\t") for i in f.readlines()]
    f.close()
    seq_id = [i[0] for i in taxa_table]
    taxonamy = [i[1] for i in taxa_table]
    pool = {}
    for i in range(len(seq_id)):
        if taxonomy_rank(rank)+taxon in taxonamy[i]:
            rep_set = otu_dict[seq_id[i]]
            current_sample = rep_set[0].split(" ")[1].split("_")[0]
            if current_sample not in pool:
                pool[current_sample] = []
            pool[current_sample].append([seq_id[i], rep_set[1]])
    samples = pool.keys()
    for sample in samples:
        f = open("origin_track/{0}.fasta".format(sample), "w")
        for each in pool[sample]:
            f.write(">" + each[0] + "\n")
            f.write(str(each[1]) + "\n")
        f.close()

def pairwise_compare():
    cmd_makeblastdb = "makeblastdb -dbtype nucl -in {0} -title {1} -out {2}"
    cmd_blast = "blastall -p blastn -d {0} -i {1} -o {2} -m 8 -b 1 -a 4"
    files = [i for i in listdir("./") if isfile(join("./",i))]
    prefixes = []
    for file in files:
        prefix = ".".join(file.split(".")[:-1])
        prefixes.append(prefix)
        os.system(cmd_makeblastdb.format(file,prefix,prefix))
    for i in range(1,len(files)):
        query = files[i]
        subjects = prefixes[:i]
        for subject in subjects:
            os.system(cmd_blast.format(subject,query,prefixes[i]+"_vs_"+subject+".tab"))

def fetch_match():
    files = [file for file in listdir("./") if isfile(join("./",file))]
    for file in files:
        print(file.split(".")[0])
        os.system("grep 100.00 {0}".format(file))
        print("\n")

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()
    working_dir = args.input
    taxon = args.taxon
    rank = args.rank
    mapping = args.mapping
    os.chdir(working_dir)
    # f = open(mapping,"r")
    # samples = [i.strip().split("\t")[0] for i in f.readlines()[1:]]
    # f.close()
    ## extract_sequences(rank=rank, taxon=taxon)
    os.chdir("origin_track/")
    ## pairwise_compare()
    os.mkdir("blast_out")
    os.system("mv *.tab blast_out/")
    os.chdir("blast_out")
    fetch_match()





# MAIN
if __name__ == "__main__":
    main()