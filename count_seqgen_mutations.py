import sys 
import re
from collections import OrderedDict

'''
Reads in a set of single nucleotides simulated on gene trees
in seq-gen with ancestral sequences written. Counts the number 
of mutations that happened on each tree.  
'''

def parse_seqgen(seqfile, ntaxa):
    '''
    Parses seq-gen file. Returns a list of ordered 
    taxon-allele pairs for each tree.   
    '''
    lines = []

    with open(seqfile) as seqs:
        for line in seqs:
            if re.match(r'\w', line): #if line starts with character
                lines.append(str.strip(line))
	
    lines = [lines[i].replace('\t', ' ') for i in range(len(lines))] #replaces tabs with space
	
    trees = [lines[i:i+(ntaxa*2)-1] for i in range(0, len(lines), (ntaxa*2)-1)] #splits into gene trees
	

    return trees	

def count_mutations(tree, ntaxa):
    '''
    Takes pairs of taxa/nodes and alleles, 
    and returns the number of mutations that 
    happened along the tree. Dictionary must be ordered
    in same way as seq-gen output. 
    '''	

    labels = [int(tree[i].split()[0]) for i in range(len(tree))] #node/taxon IDs
    alleles = [tree[i].split()[1] for i in range(len(tree))] #corresponding alleles 
    root = ntaxa + 1 #node label for root of the whole tree
    comparisons = [] #to keep track of which comparisons have already been made 
    current_taxon = 1 #to keep track of which taxon we are currently following 
    mutations = 0 #mutation counter 

    while current_taxon <= ntaxa:
        for i in range(len(labels)):

            if labels[i] > root and labels[i-1] >= root: #if the current label is an ancestral node other than the root, and the previous label is an ancestral node
                if [labels[i], labels[i-1]] not in comparisons: #if this comparison between nodes has not already been made 
                    if alleles[i] != alleles[i-1]: #if there was a mutation 
                        mutations += 1
                        comparisons.append([labels[i], labels[i-1]])
                    else:
                        comparisons.append([labels[i], labels[i-1]])

            if labels[i] == current_taxon: #if the label is the current tip taxon
                if labels[i-1] >= root: #if i - 1 is the subtending node
                    if alleles[i] != alleles[i-1]: #if there was a mutation 
                        mutations += 1
                elif labels[i-2] >= root: #if i-2 is the subtending node
                    if alleles[i] != alleles[i-2]: #if there was a mutation
                        mutations += 1
                current_taxon += 1 #update current taxon
                    
    return mutations

trees = parse_seqgen(sys.argv[1], 15)
print('Tree', '# of mutations')
for index, tree in enumerate(trees):
	print(index, count_mutations(tree, 15))
