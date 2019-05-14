import sys
from Bio import Phylo as phy
from io import StringIO

with open(sys.argv[1]) as newickfile:
    tree = newickfile.readline() #read tree as string 
    test_newick_str = tree #string version of tree
    test_newick_phy = phy.read(StringIO(tree), 'newick') #treemixin class object


def newick2splits(newicktree):
    '''
    Takes an ultrametric tree with branch lengths
    in newick format, and returns a list of ordered
    split times to pass to ms
    '''
    
    distances_pairwise = [] #for storing pairwise distances between taxa
    sisters = [] #for storing pairs of sister species
    #Get taxa list from tree string
    taxalist = [test_newick_str[i-1] for i in range(len(test_newick_str)) if test_newick_str[i] == ':' and test_newick_str[i-1].isdigit()] 
    nsplits = len(taxalist)-1 #number of taxa and number of splits 
    
    for i in range(len(taxalist)):
        for j in range(i+1, len(taxalist)): #for each pairwise comparison between taxa
            #store pairwise distances
            distances_pairwise.append([taxalist[i], taxalist[j], test_newick_phy.distance(taxalist[i], taxalist[j])])       
            
            if len(test_newick_phy.trace(taxalist[i], taxalist[j])) == 2: #if the pair are sister taxa
                sisters.append([taxalist[i], taxalist[j]]) #store in list of sister taxa
 
    splits = [] #list to store the splits to be specified to ms 

    smallers = [min(sisters[i]) for i in range(len(sisters))] #get the smaller taxon label from each sister pair 

    for i in range(len(smallers)):
        for j in range(i+1, len(smallers)): #for each pairwise comparison of ancestral populations 
            splits.append([smallers[i], smallers[j]]) #append the split 
    
    for i in range(len(sisters)): #store the splits between sister taxa
        splits.append(sisters[i])
   
      
    for i in range(len(taxalist)-1): #store the rest of the splits 
        if [sorted(taxalist)[i], sorted(taxalist)[i+1]] not in splits and [sorted(taxalist)[i+1], sorted(taxalist)[i]] not in splits:
            if len(splits) < nsplits:
                splits.append([sorted(taxalist)[i], sorted(taxalist)[i+1]])
            
    unique_pairwise = [] #to store unique and relevant pairwise distances 
    
    for i in range(len(distances_pairwise)):
        for j in range(len(splits)):
            if all(elem in distances_pairwise[i] for elem in splits[j]): #if the pairwise distance contains a relevant split 
                if distances_pairwise[i] not in unique_pairwise: #if it is a previously unseen split
                    unique_pairwise.append(distances_pairwise[i])

    split_times = sorted(unique_pairwise, key = lambda x: int(x[2]), reverse = True) #sort the list by descending order of split times
    
    for i in range(len(split_times)): #convert split times to proper units 
        split_times[i][2] = (split_times[i][2])/2
    
    for i in range(len(split_times)): #fix order of species labels
        if int(split_times[i][0]) == int(split_times[i][1]) - 1:
            split_times[i][0], split_times[i][1] = split_times[i][1], split_times[i][0]

    return split_times   

test = newick2splits(test_newick_str)
print(test)             
