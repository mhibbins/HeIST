"""
Scripts that do the parsing of gene trees and sequences.
"""

from itertools import zip_longest
from itertools import combinations
from Bio.Phylo.Consensus import _BitString
from Bio import Phylo
import io
import re
import multiprocessing as mp
import operator

resultss = 0
#disc_g = []
#conc_g = []

def grouper(iterable, n, fillvalue=None):
    """
    Utility function
    """
    args = [iter(iterable)] * n
    return(zip_longest(*args, fillvalue=fillvalue))

def cluster(d):
    """
    Utility function
    """
    clusters = {}
    for key, val in d.items():
        clusters.setdefault(val, []).append(key)
    return(clusters)

def getSisters(tree, t='g'):
    """Some nasty regex to get pairs of sister taxa (only at terminal branches)"""
    if t=='s':
        l = re.findall("\(([1-9][0-9]|\d),([1-9][0-9]|\d)\)", tree)
    else:
        l = re.findall("\(([1-9][0-9]|\d):\d\.\d\d\d,([1-9][0-9]|\d):\d\.\d\d\d\)", tree)
    return(l)
    
def checkEqual(lst):
    """
    Utility function
    """
    return(lst[1:] == lst[:-1])

def readSeqs(seqs, ntaxa, speciesPattern, nodes, batch):
    """
    Reads in sequences, determines if gene tree site pattern matches species tree
    site pattern. Returns indices of those which do.
    """
    indices = []
    c = cluster(speciesPattern)
    shouldMatch1 = c['0']
    shouldMatch2 = c['1']

    
    tmpFocal = open('focaltrees' + str(batch) + '.tmp', 'w')

    index = 0
    with open(seqs, 'rU') as f:
        for lines in grouper(f, ntaxa+nodes+1, ''):
            assert len(lines) == ntaxa+nodes+1
            pattern = {}
            index += 1
            for x, line in enumerate(lines):
                if x != 0:
                    l = line.replace('\n', '').split()
                    pattern[l[0]] = l[1]
            levels = set()
            for key, val in pattern.items():
                if int(key) in range(1,ntaxa+1):
                    levels.add(val)
            if len(levels) == 2:
                a = []
                for taxa in shouldMatch1:
                    a.append(pattern[taxa])
                if (checkEqual(a)):
                    b = []
                    for taxa in shouldMatch2:
                        b.append(pattern[taxa])
                    if (checkEqual(b)):
                        indices.append(index)
                        for y in lines:
                            tmpFocal.write(y)
    tmpFocal.close()
    return(indices)
                    
def getTrees(treefile, matchlist):
    """
    Returns list of trees at indices obtained from readSeqs
    """
    focal_trees = []
    all_trees = []
    trees_dont_follow = []
    trees = open(treefile, 'r')
    for i, line in enumerate(trees):
        l = line.replace('\n','')
        if len(l) > 3:
            all_trees.append(l)
    trees.close()
    for i, tree in enumerate(all_trees):
        if i in matchlist:
            focal_trees.append(tree)
        else:
            trees_dont_follow.append(tree)
    return(focal_trees, trees_dont_follow)

def _bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString(''.join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs

def rev(sis):
    """Utility function"""
    return((sis[1],sis[0]))

def calcDistance(tree):
    distances = {}
    for x,y in combinations(tree.get_terminals(), 2):
        distances[(x.name, y.name)] = round(tree.distance(x,y),2)
    return(distances)

def judgeDistances(speciesDistances, focalDistances):
    """
    Returns true if distances between taxa match up. 
    """
    #TODO: These can be calculated prior to calling this function to speed up...
    sp_d = {v: k for k, v in speciesDistances.items()}
    sp = {v: k for k, v in sp_d.items()}
    
    #This is the dictionary containing unique distances. They establish the rules
    #for coales times
    sorted_spdis = sorted(sp.items(), key=operator.itemgetter(1))[::-1]
    
    focal = {}
    for key, val in focalDistances.items():
        if (key in sp.keys()) or (rev(key) in sp.keys()):
            focal[key] =val
    focal = sorted(focal.items(), key=operator.itemgetter(1))[::-1]

    print(sorted_spdis)
    print(focal)

    for i in range(0, len(sorted_spdis)):
        spp = sorted_spdis[i][0]
        gene = focal[i][0]
        print(spp, gene)
        if (spp == gene) or (rev(spp) == gene) or (spp == rev(gene)) or (rev(spp) == rev(gene)):
            continue
        else:
            return(False)
    return(True)

    #Now need to check that 

def compareToSpecies(tree1, tree2, spp_sisters, species_distances):
    """Compares tree topologies. Will first check if sister taxa in the species tree are also sister 
    in the gene tree, returning false at the first non-shared occurence. If all sister taxa are present,
    it will calculate a bitstring distance with Biopython Phylo."""
    sisters = getSisters(tree2)
    top = bool
    for s in sisters:
        if (s not in spp_sisters) and (rev(s) not in spp_sisters):
            return(False)
    
    tree1 = tree1.replace(";", '')
    tree2 = tree2.replace(';', '')
 
    tree1 = Phylo.read(io.StringIO(tree1), "newick")
    tree2 = Phylo.read(io.StringIO(tree2), "newick")
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        top = False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        top = True
    else:
        top = False

    distances = calcDistance(tree2)
    dists = bool
    #print(species_distances)
    if judgeDistances(species_distances, distances) and (top == True):
        return(True)
    else:
        return(False)


def propDiscordant(focal_trees, species_tree):
    """
    Original function
    Determines the proportion of focal_trees (which have the same site pattern as the
    species tree) which are discordant (i.e. have a different topology)
    """
    i = 0
    countDis = 0
    disc_g = []
    conc_g = []

    speciesDistances = calcDistance(Phylo.read(io.StringIO(species_tree), "newick"))

    spp_sisters = getSisters(species_tree,'g')
    for i, tree in enumerate(focal_trees):
        r = call(species_tree, tree, spp_sisters, i, speciesDistances)
        if r[0] == 1:
            disc_g.append(r[1])
            countDis += 1
        elif r[0] == 0:
            conc_g.append(r[1])
    try:
        return([countDis, len(focal_trees), countDis/len(focal_trees)], disc_g, conc_g)
    except:
        return([countDis, len(focal_trees), 0.0], disc_g, conc_g)

def call(species_tree, tree, spp_sisters, i, species_distances):
    """Function to make parallel calling easier"""
    if compareToSpecies(species_tree, tree,spp_sisters, species_distances) == False:
        return([1, i])
    else:
        return([0, i])

def propDiscordant_async(focal_trees, species_tree):
    """
    Asynchronous pooling
    Determines the proportion of focal_trees (which have the same site pattern as the
    species tree) which are discordant (i.e. have a different topology)
    """
    global resultss
    global disc_g
    global conc_g
    pool = mp.Pool(mp.cpu_count())
    countDis = 0
    spp_sisters = getSisters(species_tree,'s')
    #print(spp_sisters)
    for i, tree in enumerate(focal_trees):
        pool.apply_async(call, args=(species_tree, tree, spp_sisters, i),  callback=collect_result)
    
    pool.close()
    pool.join()

    countDis = resultss
    d = disc_g
    c = conc_g
    
    try:
        resultss = 0
        disc_g = []
        conc_g = []
        return([countDis, len(focal_trees), countDis/len(focal_trees)], d, c)
    except:
        resultss = 0
        disc_g = []
        conc_g = []
        return([countDis, len(focal_trees), 0.0], d, c)

def collect_result(result):
    """Callback function for asynchronous pooling"""
    global resultss
    global disc_g
    global conc_g
    if result[0] == 1:
        disc_g.append(result[1])
    elif result[0] == 0:
        conc_g.append(result[1])
    resultss += result[0]

def parse_seqgen(seqfile, ntaxa, mask):
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
    return [trees[i] for i in mask]

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
