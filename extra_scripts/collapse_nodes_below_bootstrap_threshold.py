#!/usr/bin/env python

"""
Input: a dir of newick trees 
Output: newick trees with with nodes collapsed that
	fall below a support threshold
"""

import os,sys
import subprocess
import argparse
from ete3 import Tree, TreeNode

# create ArgumentParser object to hold all necessary info to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument("treefile", help="input directory")
parser.add_argument("threshold", help="bootstrap threshold to collapse branches")

args = parser.parse_args()

# read in treefile
threshold = args.threshold
tree = Tree(args.treefile)

for node in tree.traverse():
	if not node.is_leaf() and not node.is_root():
		if node.support < float(threshold):
				node.delete()
tree.write(format=2, outfile=args.treefile+".collapse33")
