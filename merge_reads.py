#!/usr/bin/python
"""Sometimes we want to merge raw reads from different batches into one, this process will have the re-assignment of
    barcodes to avoid duplicates.
    Files to prepare beforehand: a mapping file indicating samples, original barcodes used, which file for each sample
    belongs to,
    NAME THESE COLUMNS: SampleID,BarcodeSequence,file (This is vital, and you can name other columns whatever you like,
    but for your convenience, name them according to the standard of mapping file of QIIME)
    Save them as *.csv format.
    ;
    raw read files
    Put them all in the same directory and run the script.

"""

# IMPORT
import argparse
import pandas as pd
from Bio import SeqIO
from os.path import join
from itertools import product
import random
import sys


# FUNCTIONS
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="Sometimes we want to merge raw reads from different batches into one, "
                    "this process will have the re-assignment of barcodes to avoid duplicates.\n "
                    "Files to prepare beforehand: a tab-delimited mapping file indicating: samples, original barcodes, absolute file path "
                    "for each sample belongs to.\n "
                    "Example header line of a mapping file: (with only the necessary columns for merging, not future microbiome analysis)\n"
                    "SampleID\tBarcodeSequence\tLocation\n"
                    "Note that DO NOT put \# in the header line of this file.\n"
    )
    parser.add_argument("-m", dest="mapping", help="The mapping file name")
    parser.add_argument("-o", dest="output", help="The merged reads file name")
    args = parser.parse_args()
    return args

def merge(mapping, output):
    mapping_df = pd.read_csv(mapping)
    raw_read_files = list(set(mapping_df["Location"]))
    raw_read_dict = {}
    for raw_read_file in raw_read_files:
        f = open(raw_read_file,"r")
        raw_read_dict[raw_read_file] = list(SeqIO.parse(f,"fasta"))
        f.close()
    num_of_samples = len(mapping_df["SampleID"])
    bases = ['A', 'T', 'C', 'G']
    kmers = [''.join(i) for i in product(bases, repeat=8)]
    new_barcodes = random.sample(kmers, num_of_samples)
    mapping_df["new_barcode"] = new_barcodes
    output_handler = open(output,"w")
    indices = mapping_df.index
    barcode_translate = {}
    for i in indices:
        if mapping_df.get_value(i,"Location") not in barcode_translate:
            barcode_translate[mapping_df.get_value(i, "Location")] = {}
            barcode_translate[mapping_df.get_value(i,"Location")][mapping_df.get_value(i,"BarcodeSequence")] = \
                mapping_df.get_value(i,"new_barcode")
        else:
            barcode_translate[mapping_df.get_value(i, "Location")][mapping_df.get_value(i, "BarcodeSequence")] = \
                mapping_df.get_value(i, "new_barcode")
    print(barcode_translate)
    for raw_read_file in raw_read_files:
        for each_read in raw_read_dict[raw_read_file]:
            seq = str(each_read.seq)
            try:
                new_seq = barcode_translate[raw_read_file][seq[:8]] + seq[8:]
                output_handler.write(">" + each_read.id + "\n")
                output_handler.write(new_seq + "\n")
            except:
                continue
    output_handler.close()
    mapping_df.to_csv("new_" + mapping)



# MAIN
if __name__ == '__main__':
    args = get_parsed_args()
    mapping = args.mapping
    output = args.output
    merge(mapping=mapping, output=output)