#!/usr/bin/env python

import dendropy
import sys
import os
import copy
import argparse
import os.path

# create ArgumentParser object to hold all info necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument("treefile", help="tree in newick format")
parser.add_argument("tips_to_remove", help="list of tips to prune from treefile")
parser.add_argument("outdir", help="output directory")

args = parser.parse_args()
	
# read in treefile and tips_to_remove
treeName = args.treefile
samples = open(args.tips_to_remove)
excluded = [s[:-1] for s in samples.readlines()]
tree = dendropy.Tree.get(path=treeName, schema='newick', preserve_underscores=True)

# prune tips and write new treefile
tree.prune_taxa_with_labels(excluded)
tree.write(path=args.outdir+'/'+treeName+'.pruned', schema='newick')
	
	
