#! /usr/bin/env python3

import argparse
import re
import sys,os
import csv

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "summarize the tips removed by TreeShrink" )

# add arguments
parser.add_argument( "ending", help="file ending for summary files" )
parser.add_argument( "directory", help="working directory" )

args = parser.parse_args()

def get_taxon(name):
	name = re.sub("CCMP495_Minidiscus_variabilis_.*?_DN", "CCMP495_Minidiscus_variabilis_DN", name, flags=re.DOTALL)
	
	rnaseq_re = re.compile("(\S+)_TRINITY_\S+$")
	jgi_re = re.compile("(\S+)_XP_\S+$")
	ecto = re.compile("Ec32_Ectocarpus_siliculosus")
	fracy = re.compile("CCMP1102_Fragilariopsis_cylindrus")
	phaeo = re.compile("CCAP1005_Phaeodactylum_tricornutum")
	mini = re.compile("CCMP495_Minidiscus_variabilis")
	oceanica = re.compile("CCMP1005_Thalassiosira_oceanica")
	cryptica = re.compile("CCMP332_Cyclotella_cryptica")
	profunda = re.compile("ECT2AJA-044_Thalassiosira_profunda")
	chaetoceros = re.compile("NIES-3715_Chaetoceros_tenuissimus")
	
	rnaseq_match = re.match(rnaseq_re, name)
	jgi_match = re.match(jgi_re, name)
	ecto_match = re.match(ecto, name)
	fracy_match = re.match(fracy, name)
	phaeo_match = re.match(phaeo, name)
	mini_match = re.match(mini, name)
	oceanica_match = re.match(oceanica, name)
	cryptica_match = re.match(cryptica, name)
	profunda_match = re.match(profunda, name)
	chaetoceros_match = re.match(chaetoceros, name)
	
	if rnaseq_match:
		return rnaseq_match.group(1)
	elif jgi_match:
		return jgi_match.group(1)
	elif ecto_match:
		return "Ec32_Ectocarpus_siliculosus"
	elif fracy_match:
		return "CCMP1102_Fragilariopsis_cylindrus"
	elif phaeo_match:
		return "CCAP1005_Phaeodactylum_tricornutum"
	elif mini_match:
		return "CCMP495_Minidiscus_variabilis"
	elif oceanica_match:
		return "CCMP1005_Thalassiosira_oceanica"
	elif cryptica_match:
		return "CCMP332_Cyclotella_cryptica"
	elif profunda_match:
		return "ECT2AJA-044_Thalassiosira_profunda"
	elif chaetoceros_match:
		return "NIES-3715_Chaetoceros_tenuissimus"
	else:
		return name.split("_WR")[0]

def get_front_names(tips):
	return [get_taxon(i) for i in tips]

outfile = open("treeshrink_summary", "w")
outfile.write("taxonID\ttips_removed\n")

DIR = args.directory + "/"
DICT = {}

for i in os.listdir(DIR):
	if i[-len(args.ending):] == args.ending in i:
		if os.path.getsize(i) > 1:
			print(i)
			with open(DIR+i, "r") as infile:
				values = infile.read().split('\t')
				print(values)
			names = list(map(get_taxon, values))
			names[-1] = names[-1].strip()
			del names[-1]
			print(names)
			for taxon in names:
				if taxon not in DICT:
					DICT[taxon] = 0
				DICT[taxon] += 1
for taxon in DICT:
	outfile.write(taxon + '\t' + str(DICT[taxon]) + '\n')
outfile.close()

