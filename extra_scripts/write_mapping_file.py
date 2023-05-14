#! /usr/bin/env python3

import argparse
import re
import sys,os
import csv
from ete3 import Tree

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "create a gene name / species name mapping file" )

# add argumeents
parser.add_argument( "ending", help="file ending for alignment files" )
parser.add_argument( "directory", help="working directory" )

args = parser.parse_args()

def get_taxon(name):
	
	rnaseq_re = re.compile("(\S+)_TRINITY_\S+$")
	jgi_re = re.compile("(\S+)_XP_\S+$")
	oceanica = re.compile("CCMP1005_Thalassiosira_oceanica")
	cryptica = re.compile("CCMP332_Cyclotella_cryptica")
	profunda = re.compile("ECT2AJA-044_Thalassiosira_profunda")
	
	rnaseq_match = re.match(rnaseq_re, name)
	jgi_match = re.match(jgi_re, name)
	oceanica_match = re.match(oceanica, name)
	cryptica_match = re.match(cryptica, name)
	profunda_match = re.match(profunda, name)
	
	if rnaseq_match:
		return rnaseq_match.group(1)
	elif jgi_match:
		return jgi_match.group(1)
	elif oceanica_match:
		return "CCMP1005_Thalassiosira_oceanica"
	elif cryptica_match:
		return "CCMP332_Cyclotella_cryptica"
	elif profunda_match:
		return "ECT2AJA-044_Thalassiosira_profunda"
	else:
		return name.split("_WR")[0]

DIR = args.directory + "/"

for i in os.listdir(DIR):
	if i[-len(args.ending):] == args.ending in i:
		if os.path.getsize(i) > 1:
			print(i)
			outfile = open(i + '.mapping.txt', 'w')
			t = Tree(i)
			leaves = t.get_leaves()

			for l in leaves:
				species = get_taxon(l.name)
				outfile.write(l.name + '\t' + species + '\n')
			outfile.close()	
			
