#! /usr/bin/env python3

# this script opens BLAST output, which has been filtered ahead of time to include only the query sequences that need to be reverse complemented
# a FASTA file, potentially containing a mix of properly and improperly oriented sequences, is opened, and reverse complements output for
# all the queries in the BLAST output file


# load required modules
import sys
import os
import re
import argparse
import csv
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument( "blast", help="BLAST outfmt 6 with queries to be reverse complemented" )

args = parser.parse_args()

# keep track of number of identical query-sbjct matches
identical = 0

# keep track of number of query-sbjct matches that are the same length
same_length = 0

# keep track of number of query-sbjct matches that are full-length alignments
full_length = 0

# keep track of number of perfect matches
perfect = 0

# keep track of number of imperfect matches
imperfect = 0

# same length with similarity <100 but >90
same_length_high_sim = 0


# open and parse BLAST output
# BLAST format: "6 qseqid sseqid pident qlen slen qcovs qstart qend sstart send evalue bitscore"
with open(args.blast, 'r') as blast_file:
	for line in csv.reader(blast_file, delimiter = '\t' ):
        # are query and sbjct the same length?
		if int(line[3]) == int(line[4]):
			same_length += 1
        # are query and sbjct 100% identical?
		if float(line[2]) == 100:
			identical += 1
       	# is this a full-length alignment?
		if int(line[5]) == 100:
			full_length += 1
		# check if equal length and high similarity
		if int(line[3]) == int(line[4]) and float(line[2]) >= 90:
			same_length_high_sim += 1
       	
       	# check if it matches all three criteria
		if int(line[3]) == int(line[4]) and float(line[2]) == 100 and int(line[5]) == 100:
			perfect   += 1
		else:
			imperfect += 1
blast_file.close()

print("perfect matches:", perfect)
print("imperfect:", imperfect)
print("100% identity:", identical)
print("same-length:", same_length)
print("full-length matches:", full_length)
print("same-len-high-sim matches:", same_length_high_sim)
