#! /usr/bin/python3

import argparse
import sys,os
import math
import re
import copy
import numpy as np
import numpy.polynomial.polynomial as poly
from ete3 import Tree

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "Transform branch lengths to coalescent units" )

# add arguments
parser.add_argument( "--tree", "-t", help="Tree file annotated with quartet concordance factor)" )

args = parser.parse_args()

def main():
	newick = open(args.tree, 'r').read()
	
	treeSp,t,treeSp_lower_CI,t_lower_CI,treeSp_upper_CI,t_upper_CI,intercept,coef,n,c = subs2coal(newick)
	print(treeSp)

# begin calculations
def subs2coal(newick_string):
        
        '''
        Takes a newick string with the nodes labelled with concordance factors,
        and returns the same string with branch lengths converted to coalescent
        units.
        '''

        cfs = re.findall("\)(.*?)\:", newick_string) #regex to get concordance factors
        cfs = list(filter(None, cfs)) #removes empty values (root has no scf)
        coal_internals = []

        for i in range(len(cfs)):

                if (float(cfs[i])/100) < 1 and float(cfs[i]) != 1: #catch for missing data / values of 1
                        coal_estimate = -1*(math.log(3/2) + math.log(1 - float(cfs[i])/100))
                        if coal_estimate > 0:
                                coal_internals.append(coal_estimate)
                        else:
                                coal_internals.append(0.01)
                else:
                        coal_internals.append(np.NaN)

        coal_internals = [float(i) for i in coal_internals]
        newick_internals = []
        internal_sections = []
        sections = newick_string.split(':')

        for i in range(len(sections)):
                if len(sections[i].split(")")) > 1 and sections[i].split(")")[1] in cfs:
                                newick_internals.append(newick_string.split(':')[i+1].split(',')[0].split(")")[0])
                                internal_sections.append(sections[i+1])

        
        newick_internals = [float(i) for i in newick_internals]
         
        def branch_regression(newick_branches, coal_internals):

                '''
                Returns the regression coefficient of the coalescent branch lengths
                on the branch lengths in the newick string
                '''

                newick_reg_branches = copy.deepcopy(newick_branches)
                coal_reg_internals = copy.deepcopy(coal_internals)

                for i in range(len(newick_branches)): #drops missing data
                        if np.isnan(newick_branches[i]) or np.isnan(coal_internals[i]):
                                newick_reg_branches[i] = "N/A"
                                coal_reg_internals[i] = "N/A"
                                #del newick_branches[i]
                                #del coal_internals[i]

                while "N/A" in newick_reg_branches:
                    newick_reg_branches.remove("N/A")

                while "N/A" in coal_reg_internals:
                    coal_reg_internals.remove("N/A")

                #print(newick_reg_branches)
                #print(coal_reg_internals)

                intercept, slope = poly.polyfit(newick_reg_branches, coal_reg_internals, 1)

                return intercept, slope

        intercept, coef = branch_regression(newick_internals, coal_internals)
        n = newick_internals
        c = coal_internals
        tip_sections = []

        for i in range(len(sections)): #gets the sections with tip lengths
                if sections[i] not in internal_sections:
                        tip_sections.append(sections[i])

        newick_tips = []

        for i in range(len(tip_sections)): #gets the tip lengths
                if len(tip_sections[i].split(',')) > 1:
                        newick_tips.append(tip_sections[i].split(',')[0])
                elif len(tip_sections[i].split(')')) > 1:
                        newick_tips.append(tip_sections[i].split(')')[0])
        
        newick_tips = [float(i) for i in newick_tips]

        coal_tips = [(newick_tips[i]*coef + intercept) for i in range(len(newick_tips))]

        for i in range(len(coal_tips)):
                if coal_tips[i] <= 0:
                        coal_tips[i] = 0.01
                else:
                        coal_tips[i] = coal_tips[i]

        prediction_stdev = np.std(coal_tips) #standard deviation of tip predictions
        coal_tips_lower_CI = [(tip - 1.96*prediction_stdev) for tip in coal_tips]
        coal_tips_upper_CI = [(tip + 1.96*prediction_stdev) for tip in coal_tips]

        coal_internals = [(newick_internals[i]*coef + intercept) if np.isnan(coal_internals[i]) else coal_internals[i] for i in range(len(newick_internals))]
        coal_internals_lower_CI = [((newick_internals[i]*coef + intercept) - 1.96*prediction_stdev) if np.isnan(coal_internals[i]) else coal_internals[i] for i in range(len(newick_internals))]
        coal_internals_upper_CI = [((newick_internals[i]*coef + intercept) + 1.96*prediction_stdev) if np.isnan(coal_internals[i]) else coal_internals[i] for i in range(len(newick_internals))]

        lengths = re.findall("\d+\.\d+", newick_string)

        lengths = [lengths[i] for i in range(len(lengths)) if float(lengths[i]) < 1]

        coal_lengths = []
        coal_lengths_lower_CI = []
        coal_lengths_upper_CI = []

        for i in range(len(lengths)):
                if newick_internals.count(lengths[i]) > 1 or newick_tips.count(lengths[i]) > 1:
                        sys.exit('Error: Duplicate branch lengths')
                elif float(lengths[i]) in newick_internals:
                        coal_lengths.append(coal_internals[newick_internals.index(float(lengths[i]))])
                        coal_lengths_lower_CI.append(coal_internals_lower_CI[newick_internals.index(float(lengths[i]))])
                        coal_lengths_upper_CI.append(coal_internals_upper_CI[newick_internals.index(float(lengths[i]))])
                elif float(lengths[i]) in newick_tips:
                        coal_lengths.append(coal_tips[newick_tips.index(float(lengths[i]))])
                        coal_lengths_lower_CI.append(coal_tips_lower_CI[newick_tips.index(float(lengths[i]))])
                        coal_lengths_upper_CI.append(coal_tips_upper_CI[newick_tips.index(float(lengths[i]))])
                elif float(lengths[i]) == 0: #deals with roots of length 0 in smoothed trees
                        coal_lengths.append(float(0))
        
        coal_lengths = [str(coal_lengths[i]) for i in range(len(coal_lengths))]
        coal_lengths_lower_CI = [str(coal_lengths_lower_CI[i]) for i in range(len(coal_lengths_lower_CI))]
        coal_lengths_upper_CI = [str(coal_lengths_upper_CI[i]) for i in range(len(coal_lengths_upper_CI))]
       
        coal_newick_string = copy.deepcopy(newick_string)
        coal_newick_string_lower_CI = copy.deepcopy(newick_string)
        coal_newick_string_upper_CI = copy.deepcopy(newick_string)

        for i in range(len(lengths)):
             coal_newick_string = coal_newick_string.replace(str(lengths[i]), str(coal_lengths[i]))
             coal_newick_string_lower_CI = coal_newick_string_lower_CI.replace(str(lengths[i]), str(coal_lengths_lower_CI[i]))
             coal_newick_string_upper_CI = coal_newick_string_upper_CI.replace(str(lengths[i]), str(coal_lengths_upper_CI[i]))
    
        cfs = [float(x) for x in cfs]
        cfs2 = []
        for val in cfs:
            if val.is_integer():
                cfs2.append(str(val) + ".0")
            else:
                cfs2.append(str(val))


        for i in range(len(cfs)):
                coal_newick_string = coal_newick_string.replace(str(cfs[i]), '')
                coal_newick_string_lower_CI = coal_newick_string_lower_CI.replace(str(cfs[i]), '')
                coal_newick_string_upper_CI = coal_newick_string_upper_CI.replace(str(cfs[i]), '')

        return(coal_newick_string, Tree(coal_newick_string, format=1), coal_newick_string_lower_CI, Tree(coal_newick_string_lower_CI, format=1),
                coal_newick_string_upper_CI, Tree(coal_newick_string_upper_CI, format=1), intercept, coef, n, c)

# call main function
if __name__ == '__main__':
	main()

