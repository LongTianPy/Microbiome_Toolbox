#!/usr/bin/python
"""
"""

# IMPORT
from Bio import SeqIO
import sys

# FUNCTIONS
def parse_blast_out(gene):
    blast_out = gene + ".txt"
    f = open(gene+".fas","r")
    seqs = list(SeqIO.parse(f,"fasta"))
    f.close()
    length = len(seqs[0].seq)
    f = open(blast_out,"r")
    lines = [i.strip().split("\t") for i in f.readlines()]
    f.close()
    alignment = {}
    for line in lines:
        subject = line[1]
        start = int(line[8])
        end = int(line[9])
        if subject not in alignment and int(subject) <= 410 and int(line[3])==length:
            alignment[subject] = [start,end]
        else:
            continue
    return alignment

def get_seq(alignment,db):
    f = open(db,"r")
    seqs = list(SeqIO.parse(f,"fasta"))
    f.close()
    seqs_dict = {seq.id:seq.seq for seq in seqs}
    genome_ids = alignment.keys()
    genes = []
    for id in genome_ids:
        start = alignment[id][0]
        end = alignment[id][1]
        if end > start:
            seq = str(seqs_dict[id][start-1:end])
        else:
            seq = str(seqs_dict[id][end-1:start].reverse_complement())
        genes.append([id,seq])
    return genes

def fetch_seq(gene,db):
    alignment = parse_blast_out(blast_out=blast_out)
    genes = get_seq(alignment=alignment,db=db)
    f = open(gene+"_align.fasta","r")
    for i in genes:
        f.write(">"+i[0]+"\n")
        f.write(i[1]+"\n")
    f.close()

# MAIN
if __name__ == '__main__':
    gene = sys.argv[1]
    db = sys.argv[2]
    fetch_seq(gene,db)
