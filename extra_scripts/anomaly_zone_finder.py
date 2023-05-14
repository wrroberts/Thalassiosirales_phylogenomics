#! /usr/bin/python3

import sys,os
import math
from math import exp
import dendropy
from dendropy.utility import bitprocessing
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "Test branches for the anomaly zone" )

# add arguments
parser.add_argument( "--tree", "-t", help="Tree file (branches lengths in coalescent units)" )
parser.add_argument( "--bootstraps", "-b", help="Tree file containing bootstrap trees (optional)" )
parser.add_argument( "--format", "-f", help="Format of tree file: 'nexus' or 'newick' (default: newick)", default="newick" )
parser.add_argument( "--output", "-o", help="Output file name for annotation tree (default: out.tre)", default="out.tre" )

args = parser.parse_args()

# begin calculations
def main():
	params = args
	print("\n##### anomaly_finder.py #####\n")
	# read trees
	taxon_set = dendropy.TaxonNamespace()
	print("Reading species tree...", args.tree)
	mrc = dendropy.Tree.get_from_path(src=args.tree, schema=args.format, taxon_namespace=taxon_set)
	if args.bootstraps:
		print("Reading bootstrap trees...", args.bootstraps)
		trees = dendropy.TreeList.get_from_path(args.bootstraps, args.format, taxon_namespace=taxon_set)
	else:
		print("No bootstrap trees provided.")
		
	# initialize dictionaries
	master_dict = dict()
	mrc_dict = dict()
	split_dict = dict()
	split_occ_dict = dict()
	
	# if bootstrap trees provided, get results for them
	# master_dict:
	# key = pair of edges
	# value = list of anomaly zone results (1=True; 0=False)
	if args.bootstraps:
		print("Finding anomalous nodes in bootstrap trees...")
		for sp_tree in trees:
			sp_tree.encode_bipartitions()
			master_dict = get_nodes(sp_tree, master_dict)
			split_occ_dict = split_freq(sp_tree, taxon_set, split_occ_dict)
			
	# get results for main tree
	print("Finding anomalous nodes in species tree...")
	mrc.encode_bipartitions()
	mrc_dict = get_nodes(mrc, mrc_dict) # get results for species tree
	split_dict = split_mapper(mrc, taxon_set) # match list of taxa to split pattern from the bitmask for the whole tree
	
	print("\n##### ENCODING #####\n")
	
	print("Writing results with splits encoded as:")
	labelTree = getBitStringLabels(mrc, len(taxon_set))
	labelTree.print_plot(show_internal_node_labels=True)
	
	print("(If tree is large, import the following into FigTree)\n")
	print(labelTree.as_string(schema="newick"))
	
	# if bootstraps provided:
	# for each internal edge pair:
	print("\n##### RESULTS #####\n")
	print("Showing anomalous nodes only:\n")
	if args.bootstraps:
		print("ParentEdge\tDescendantEdge\tProp.Anomalous")
		seen = False
		for k,v in master_dict.items():
			taxlab1 = list()
			taxlab2 = list()
			if k in mrc_dict.keys():
				parentLabel = bitprocessing.int_as_bitstring(k[0], length=len(taxon_set))
				parentList = taxon_set.bitmask_taxa_list(k[0])
				
				descendantLabel = bitprocessing.int_as_bitstring(k[1], length=len(taxon_set))
				descendantList = taxon_set.bitmask_taxa_list(k[1])
				
				for i in parentList:
					taxlab1.append(i.label)
				for j in descendantList:
					taxlab2.append(j.label)
				if sum(v) > 0:
					seen = True
					print(parentLabel, descendantLabel, sum(v)/len(v), sep="\t")
					
		if seen==False:
			print("No anomalous nodes detected.")
	
	else:
		print("ParentEdge\tDescendantEdge")
		seen=False
		for k,v in mrc_dict.items():
			taxlab1 = list()
			taxlab2 = list()
			anom = False
			
			parentLabel = bitprocessing.int_as_bitstring(k[0], length=len(taxon_set))
			parentList = taxon_set.bitmask_taxa_list(k[0])
			
			descendantLabel = bitprocessing.int_as_bitstring(k[1], length=len(taxon_set))
			descendantList = taxon_set.bitmask_taxa_list(k[1])
			
			for i in parentList:
				taxlab1.append(i.label)
			for j in descendantList:
				taxlab2.append(j.label)
			if sum(v) > 0:
				seen=True
				print(parentLabel, descendantLabel, sep="\t")
				
		if seen==False:
			print("No anomalous nodes detected.")
	
	print("\n##### OUTPUT #####\n")
	
	if args.bootstraps:
		print("Tree with anomalous edge pairs annotated as \"AZ Event : Prop. BS\"")
		outTreeBS = getAnomalyLabelsBS(mrc, master_dict)
		outTreeBS.print_plot(show_internal_node_labels=True)
		
		print("(If tree is large, import the following into FigTree)\n")
		print(outTreeBS.as_string(schema="newick"))
		with open(args.output, 'w') as fh:
			l = outTreeBS.as_string(schema="newick")
			l = l+"\n"
			fh.write(l)
		fh.close()
	else:
		print("Tree with anomalous edge pairs annotated (AZ pairs numbered)")
		outTree = getAnomalyLabels(mrc, mrc_dict)
		outTree.print_plot(show_internal_node_labels=True)
		
		print("(If tree is large, import the following into FigTree)\n")
		print(outTree.as_string(schema="newick"))
		with open(args.output, 'w') as fh:
			l = outTree.as_string(schema="newick")
			l = l+"\n"
			fh.write(l)
		fh.close()
		
	print("\n##### DONE #####\n\n")
	
# function to return tree with AZ edge pairs labelled
def getAnomalyLabels(tree, d):
	tree2 = tree.clone(depth=1)
	labelDict = dict()
	az_event = 1
	for k,v in d.items():
		#if AZ test positive
		if sum(v) > 0:
			for edge in k:
				if str(edge) in labelDict:
					#print("exists")
					labelDict[str(edge)] = str(labelDict[str(edge)]) + "/" + str(az_event)
				else:
					labelDict[str(edge)] = str(az_event)
			#increment event number
			az_event += 1
	for node in tree2.levelorder_node_iter():
		if node.parent_node:
			lab=str(node.edge.split_bitmask)
			if lab in labelDict:
				l = "\"" + labelDict[str(lab)] + "\""
				node.label=l
			else:
				lab = "\"" + "0" + "\""
				node.label=lab
		else:
			lab = "\"" + "0" + "\""
			node.label=lab
	return(tree2)
	
# function to return tree with AZ edge pairs labeled, for bootstrap trees
def getAnomalyLabelsBS(tree, d):
	tree2 = tree.clone(depth=1)
	labelDict = dict()
	bsDict = dict()
	az_event = 1
	for k,v in d.items():
		#if AZ test positive
		if sum(v) > 0:
			for edge in k:
				if str(edge) in labelDict:
					#print("exists")
					labelDict[str(edge)] = str(labelDict[str(edge)]) + "/" + str(az_event)
					bsDict[str(edge)] = str(bsDict[str(edge)]) + "/" + str(sum(v)/len(v))
				else:
					labelDict[str(edge)] = str(az_event)
					bsDict[str(edge)] = str(sum(v)/len(v))
			#increment event number
			az_event += 1
	for node in tree2.levelorder_node_iter():
		if node.parent_node:
			lab=str(node.edge.split_bitmask)
			if lab in labelDict:
				node.label=labelDict[str(lab)] + ":" + bsDict[str(lab)]
			else:
				node.label="0:0"
		else:
			node.label="0:0"
	return(tree2)
	
# function to return tree with bitstring node labels
def getBitStringLabels(tree, l):
	tree2 = tree.clone(depth=1)
	for node in tree2.postorder_node_iter():
		lab = "\"" + str(bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=l)) + "\""
		node.label=lab
	return tree2
	
# for a pair of edges, return 1 if anomalous, return 0 if not anomalous
def anomaly_calc(nodes):
	lambda1=nodes[0]
	lambda2=nodes[1]
	if lambda1 is not None:
		Z1=math.log(2.0/3+((3*exp(2*lambda1)-2)/(18*(exp(3*lambda1)-exp(2*lambda1))))) #function a(x) from Degnan and Rosenberg (2006)
		if lambda2 <= Z1:
			return 1
		else:
			return 0
	else:
		return 0

def split_mapper(tree, taxon_set):
	split_dict={}
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			taxlist=[]
			tax = taxon_set.bitmask_taxa_list(node.edge.split_bitmask)
			for i in tax:
				taxlist.append(i.label)
			split_dict[bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=len(taxon_set))] = taxlist
	return split_dict

# function to split
def split_freq(tree,taxon_set, split_occ_dict):
	#for each node 
	for node in tree.postorder_node_iter():
		#if not root
		if node.parent_node is None:
			node.value = 1.0
		else:
			#split=dendropy.treesplit.split_as_string(node.edge.split_bitmask, width=len(taxon_set)) #deprecated
			split=bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=len(taxon_set))
			print(split)
			if split not in split_occ_dict.keys():
				split_occ_dict[split] = [1]
			else:
				split_occ_dict[split].append(1)
	return split_occ_dict
	
# function to test for internal edge pairs under anomaly zone
# modifies a dictionary where key=pair of edges, value=list of results (1=true; 0=false)
def get_nodes(tree, master_dict):
	pair_list=[]
	#for each internal node
	#calculate
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			#root
			node.value = 1.0
		else:
			#if not root, and internal node
			if node.taxon is None:
				#get edge length of parent, and own edge length
				#x = parent; y = current
				node_pair = (node.parent_node.edge_length, node.edge_length)
				#calculate if pair of edges are under anomaly boundary
				anomalous=anomaly_calc(node_pair)
				#print(anomalous)
				edge_pair = (node.parent_node.edge.split_bitmask, node.edge.split_bitmask)
				#for edge pair, keep anomaly zone tests results
				if edge_pair not in master_dict.keys():
					master_dict[edge_pair] = [anomalous]
				else:
					master_dict[edge_pair].append(anomalous)	
	return master_dict		

# call main function
if __name__ == '__main__':
	main()
