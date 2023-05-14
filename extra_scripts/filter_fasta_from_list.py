#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument( "fasta", help="name of FASTA file" )
parser.add_argument( "remove_list", help="name of file with list of sequences to remove" )

args = parser.parse_args()

delete_seqs = []

# open and parse file with list of sequences to exclude
exclude_file = open( args.remove_list, 'r' )

for line in exclude_file:
    delete_seqs.append(line.rstrip())

# close list file
exclude_file.close()


# open and parse file1 (FASTA file that we're filtering contaminants from)
fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')

# open output file
for record in fasta_sequences:
    if record.description not in delete_seqs:
    	# print(record.id)
    	print(record.format("fasta"))
