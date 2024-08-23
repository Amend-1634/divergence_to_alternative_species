#!/usr/bin/env python3

from sys import argv
import operator
import gzip
from itertools import islice
from Bio import SeqIO

#read fasta of reference in dictionary, allow gzipped file
if argv[1].endswith(".gz"):
    with gzip.open(argv[1], 'rt') as f1:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(f1, "fasta"))
 
else:
    with open(argv[1]) as f1:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(f1, "fasta"))

#read fasta of sample in dictionary, allow gzipped file
if argv[2].endswith(".gz"):
    with gzip.open(argv[2], 'rt') as f1:
        sample_dict = SeqIO.to_dict(SeqIO.parse(f1, "fasta"))

else:
    with open(argv[2]) as f1:
        sample_dict = SeqIO.to_dict(SeqIO.parse(f1, "fasta"))


transitions_type = ["AG","GA","CT","TC"]

#initiate transition and transversion counts
mutations = 0
transitions = 0
transversions = 0
total_sites_covered = 0

# Get the intersection of keys
intersection_keys = set(fasta_dict.keys()) & set(sample_dict.keys())

# Extract values from dict1 based on intersection keys
fasta_dict_ext = {key: fasta_dict[key] for key in intersection_keys}


#compare to get the number of transitions and transversions
for contig,seq in fasta_dict_ext.items():
    reference_sequence = seq
    sample_sequence = sample_dict[contig]
    for i in range(len(reference_sequence)):
        #get rid of the first and last 3 base pair (the end of the paired-end reads):
        if i > 2 and i < len(reference_sequence) - 3:
            #counts for mutations and total sites covered
            if sample_sequence[i].upper() != "N" and reference_sequence[i].upper() != "N":
                total_sites_covered += 1
                if reference_sequence[i].upper() != sample_sequence[i].upper():
                    mutations += 1
                    #transition and transversion
                    if reference_sequence[i].upper() + sample_sequence[i].upper() in transitions_type:
                        transitions += 1
                    else:
                        transversions += 1

print(mutations)

print(transitions)

print(transversions)

print(total_sites_covered)
