import sys 
import re

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
    happened along the tree. Pairs must be ordered
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

            if labels[i] > root and labels[i-1] >= root and labels[i] == labels[i-1] + 1: #if the current label is an ancestral node other than the root, and the previous label is an ancestral node
                if [labels[i], labels[i-1]] not in comparisons: #if this comparison between nodes has not already been made 
                    if alleles[i] != alleles[i-1]: #if there was a mutation 
                        mutations += 1
                        comparisons.append([labels[i], labels[i-1]])
                    else:
                        comparisons.append([labels[i], labels[i-1]])
            elif labels[i] > root and labels[i-2] >= root and labels[i] == labels[i-2] + 1: #if the current label is an ancestral node, and the previous ancestral node has 1 taxon descending from it
                if [labels[i], labels[i-2]] not in comparisons: #if this comparison between nodes has not already been made 
                    if alleles[i] != alleles[i-2]: #if there was a mutation 
                        mutations += 1
                        comparisons.append([labels[i], labels[i-2]])
                    else:
                      comparisons.append([labels[i], labels[i-2]])
            elif labels[i] > root and labels[i-3] >= root and labels[i] == labels[i-3] + 1: #if the current label is an ancestral node, and the previous ancestral node has 2 taxa descending from it
                if [labels[i], labels[i-3]] not in comparisons: #if this comparison between nodes has not already been made 
                    if alleles[i] != alleles[i-3]: #if there was a mutation 
                        mutations += 1
                        comparisons.append([labels[i], labels[i-3]])
                    else:
                      comparisons.append([labels[i], labels[i-3]])
            elif labels[i] > root and labels[i-4] >= root and labels[i] == labels[i-4] + 2: #if the current label is an ancestral node, and a clade with a non-subtending ancestral node is listed before it
                if [labels[i], labels[i-4]] not in comparisons: #if this comparison between nodes has not already been made 
                    if alleles[i] != alleles[i-4]: #if there was a mutation 
                        mutations += 1
                        comparisons.append([labels[i], labels[i-4]])
                    else:
                      comparisons.append([labels[i], labels[i-4]])

            elif labels[i] == current_taxon: #if the label is the current tip taxon
                if labels[i-1] >= root: #if i - 1 is the subtending node
                    if alleles[i] != alleles[i-1]: #if there was a mutation 
                        mutations += 1
                elif labels[i-2] >= root: #if i-2 is the subtending node
                    if alleles[i] != alleles[i-2]: #if there was a mutation
                        mutations += 1
                current_taxon += 1 #update current taxon

    return mutations

def get_interesting(trees, nderived, ntaxa):
    '''
    This function uses count_mutations to 
    pull out all the "interesting" cases of 
    incongruence. For now, these are cases 
    where the number of mutations is greater
    than 1 but less than the number of derived 
    taxa. 
    '''
    interesting = [] #list of interesting trees
    
    for index, tree in enumerate(trees):
        
        if count_mutations(tree, ntaxa) > 1 and count_mutations(tree, ntaxa) < nderived: 
            interesting.append(tree)
            

    return interesting 

def summarize_interesting(tree, ntaxa):
    '''
    Summarizes the mutations that have occurred
    on the given tree.
    '''
    labels = [int(tree[i].split()[0]) for i in range(len(tree))] #node/taxon IDs
    alleles = [tree[i].split()[1] for i in range(len(tree))] #alleles
    root = ntaxa + 1
    current_taxon = 1 #current taxon tracker 
    ancestral_allele = alleles[labels.index(root)] #ancestral allele for the tree
    summary = [] #list of how each taxon got its mutation 

    while current_taxon <= ntaxa:
        for i in range(len(labels)):
   
            if labels[i] == current_taxon:
                
                if alleles[i] == ancestral_allele: #if the taxon has ancestral state
                    current_taxon += 1 #move on 

                elif alleles[i] != ancestral_allele: #if taxon has the derived state
                    if labels[i-1] >= root: #if i-1 is the subtending node 
                        if alleles[i] != alleles[i-1]: #if there was a mutation on the tip branch
                            summary.append('Taxon ' + str(labels[i]) + ' mutated to the derived state')
                        else: #if the derived state was inherited from an ancestor
                            summary.append('Taxon ' + str(labels[i]) + ' inherited the derived state from ancestral node ' + str(labels[i-1]))
                    elif labels[i-2] >= root: #if i-2 is the subtending node 
                        if alleles[i] != alleles[i-2]: #if there was a tip mutation
                            summary.append('Taxon ' + str(labels[i]) + ' mutated to the derived state')
                        else: #if derived state was inherited
                            summary.append('Taxon ' + str(labels[i]) + ' inherited the derived state from ancestral node ' + str(labels[i-2]))
                    current_taxon += 1 #update current taxon  
  
    return summary

   
test_trees = parse_seqgen(sys.argv[1], 6)
test_interesting = get_interesting(test_trees, 3, 6)
test_summarize = summarize_interesting(test_interesting[1], 6)
print(test_summarize) 

