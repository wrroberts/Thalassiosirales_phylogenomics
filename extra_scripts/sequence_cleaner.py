#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

# create ArgumentParser object to hold all necessary info to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument('fasta', help='the protein fasta file to clean')
parser.add_argument('min_length', help='the minimum protein length to keep')
parser.add_argument('seq_type', help='sequence type, protein or codon', choices=['protein','codon'])

args = parser.parse_args()

# define the function
def sequence_cleaner(fasta, min_length, seq_type):
	# create a hash table to add the sequences
	sequences = {}
	
	# use biopython to parse the fasta records
	for seq_record in SeqIO.parse(fasta, 'fasta'):

		# take the current sequence
		sequence = str(seq_record.seq).upper()
		
		# check if the current sequence is okay according to user parameters
		if seq_type == 'protein':
			if ( sequence.startswith('M') and len(sequence) >= int(min_length) ):
				if sequence not in sequences:
					sequences[sequence] = seq_record.id
				else:
					continue
		else:
			if ( sequence.startswith('ATG') and len(sequence) >= int(min_length) ):
				if sequence not in sequences:
					sequences[sequence] = seq_record.id
				else:
					continue
		
	# write the clean sequences
	# create a file in the same directory
	with open(args.fasta + '.clean','w') as output_file:
		# just read the hash table and write the file as fasta format
		for sequence in sequences:
			output_file.write('>' + sequences[sequence] + '\n' + sequence + '\n')

with open(args.fasta, 'r') as fasta:
	sequence_cleaner(fasta, args.min_length, args.seq_type)
	file_name = os.path.basename(args.fasta)
	print(file_name, 'CLEAN!')
