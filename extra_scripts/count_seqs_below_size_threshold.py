#!/usr/bin/env python3

import os
import sys
import pandas as pd
import argparse

# create ArgumentParser object to hold all necessary info to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument("dir", help="directory of input")
parser.add_argument("output", help="name of output file")
parser.add_argument("percent", help="count sequences that are shorter than this percentage")

args = parser.parse_args()

DIR = args.dir+'/'
outname = args.output
outfile = open(outname, 'w')
outfile.write('orthogroupID\ttotal_seqs\tseqs_to_remove\n')
count = 0
for i in os.listdir(DIR):
	if i.endswith('.lengths'):
		orthogroupID = i.split("_alignment")[0]
		print(orthogroupID)
		count += 1
		df = pd.read_csv(i, sep='\t', header=None)
		total_seqs = df.shape[0]
		seqs_to_remove = df[df[2] <= float(args.percent)].shape[0]
		#frequencies = df[2].value_counts()
		#seqs_to_remove = (frequencies < float(args.percent)).sum()

		outfile.write(orthogroupID+'\t'+str(total_seqs)+'\t'+str(seqs_to_remove))
		outfile.write('\n')
