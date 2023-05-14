#!/usr/bin/env python3

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd
import sys
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "BLAST search a set of protein fasta files")

# add arguments
parser.add_argument( "list", help="list of protein fasta files to BLAST search" )
parser.add_argument( "-d", "--directory", help="location of the fasta files to BLAST" )
parser.add_argument( "-r", "--reference", help="reference database fasta to BLAST against" )

args = parser.parse_args()

df = pd.read_csv(args.list, sep="\t")
out = pd.read_csv(args.list, sep="\t", usecols=["Fasta"] )

out["Best_hit"]	= "NO_HIT"
out["E_val"]	= "NO_HIT"
out["iden"]		= "NO_HIT"
out["len"]		= "NO_HIT"
out["Coords"]	= "NO_HIT"

# go through each row of the file list
for row in range(0, len(df)):
	
	# fasta in question
	OG = df.iloc[row,0]	
	#print(OG)
	
	# prepare loops
	hit = False
	
	while hit == False:
	
		max_len = 0
	
		infasta = f"{args.directory}/" + OG
		print(infasta)
	
		# search only the longest sequence from each fasta
		records = list(SeqIO.parse(infasta, "fasta"))
		records.sort(key=lambda r: -len(r))

		longest = records[0]
		print(longest.id)

		# extract sequence
		SeqIO.write(longest, "tmp.fasta", "fasta")
			
		# blast sequence against reference db
		blastx_cline = NcbiblastpCommandline(query="tmp.fasta", db=args.reference, evalue=1e-6,  outfmt=5, out="my_blast.xml", num_alignments=1, soft_masking="true", lcase_masking="true", max_hsps=1, num_threads=8, seg="yes")
			
		stdout, stderr = blastx_cline()
			
		if len(stderr) > 0: print("BLAST ERROR: ", stderr)
			
		for record in NCBIXML.parse(open("my_blast.xml")):
			if record.alignments: # skip queries with no matches
				
				# record.alignments is list of alignments ordered by e-val. I only want the first element now
				align = record.alignments[0]
				hsps = align.hsps[0]
				e_val = hsps.expect
				
				iden = hsps.identities
				leng = hsps.align_length
				
				out.loc[row, "Best_hit"] = align.title
				out.loc[row, "E_val"] = e_val
				out.loc[row, "iden"] = iden
				out.loc[row, "len"] = leng
				out.loc[row, "Coords"] = str( [hsps.sbjct_start, hsps.sbjct_end] )
				
				hit = True # End while loop
				
		else: hit = True
		
# output will have same name as input list plus "_blasted.tsv" suffix
outname = args.list[:-4] + "_blasted.tsv"
out.to_csv(outname, sep="\t", index=False)
	
