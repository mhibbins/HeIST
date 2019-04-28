"""
Scripts that do the parsing of gene trees and sequences.
"""

from ete3 import Tree
from itertools import zip_longest
from Bio.Phylo.Consensus import _BitString
from Bio import Phylo
import io
import re
import multiprocessing as mp

resultss = 0

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
    """Some nasty regex"""
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

def readSeqs(seqs, ntaxa, speciesPattern):
    """
    Reads in sequences, determines if gene tree site pattern matches species tree
    site pattern. Returns indices of those which do.
    TODO: Make able to handle seq-gen output with ancestral sequences. Currently will fail.
    """
    indices = []
    c = cluster(speciesPattern)
    shouldMatch1 = c['0']
    shouldMatch2 = c['1']

    print(shouldMatch1)
    print(shouldMatch2)

    index = 0
    with open(seqs, 'rU') as f:
        for lines in grouper(f, ntaxa+1, ''):
            assert len(lines) == ntaxa+1
            pattern = {}
            index += 1
            for x, line in enumerate(lines):
                if x != 0:
                    l = line.replace('\n', '').split()
                    pattern[l[0]] = l[1]
            #print(pattern)
            levels = set()
            for val in pattern.values():
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
    return((sis[1],sis[0]))

def compareToSpecies(tree1, tree2, spp_sisters):
    sisters = getSisters(tree2)
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
        return(False)
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return(True)
    else:
        return(False)

"""
Depricated
def compareToSpecies(speciesTree, geneTree):

    speciesTree = Tree(speciesTree)
    geneTree = Tree(geneTree)
    print(speciesTree)
    print(geneTree)
    r = speciesTree.compare(geneTree, unrooted=False)['rf']
    print(r)
    if r == 0.0:
        return(True)
    else:
        #print(geneTree)
        return(False)
"""

def propDiscordant(focal_trees, species_tree):
    """
    Original function
    Determines the proportion of focal_trees (which have the same site pattern as the
    species tree) which are discordant (i.e. have a different topology)
    """
    i = 0
    countDis = 0
    spp_sisters = getSisters(species_tree,'s')
    for tree in focal_trees:
        if not compareToSpecies(species_tree, tree, spp_sisters):
            countDis += 1
        i += 1
    try:
        return(countDis, len(focal_trees), countDis/len(focal_trees))
    except:
        return(countDis, len(focal_trees), 0.0)

def call(species_tree, tree, spp_sisters):
    """Function to make parallel calling easier"""
    if not compareToSpecies(species_tree, tree,spp_sisters):
        return(1)
    else:
        return(0)

def propDiscordant_par(focal_trees, species_tree):
    """
    Synchronous pooling
    Determines the proportion of focal_trees (which have the same site pattern as the
    species tree) which are discordant (i.e. have a different topology)
    """

    pool = mp.Pool(mp.cpu_count())
    i = 0
    countDis = 0
    spp_sisters = getSisters(species_tree,'s')
    
    results = [pool.apply(call, args = (species_tree, tree, spp_sisters)) for tree in focal_trees]
    pool.close()

    countDis = sum(results)
    #print(results)

    try:
        return(countDis, len(focal_trees), countDis/len(focal_trees))
    except:
        return(countDis, len(focal_trees), 0.0)

def propDiscordant_async(focal_trees, species_tree):
    """
    Asynchronous pooling
    Determines the proportion of focal_trees (which have the same site pattern as the
    species tree) which are discordant (i.e. have a different topology)
    """
    global resultss
    pool = mp.Pool(mp.cpu_count())
    i = 0
    countDis = 0
    spp_sisters = getSisters(species_tree,'s')

    for tree in focal_trees:
        pool.apply_async(call, args=(species_tree, tree, spp_sisters),  callback=collect_result)
    
    pool.close()
    pool.join()

    countDis = resultss
    
    try:
        resultss = 0
        return(countDis, len(focal_trees), countDis/len(focal_trees))
    except:
        resultss = 0
        return(countDis, len(focal_trees), 0.0)

def collect_result(result):
    """Callback function for asynchronous pooling"""
    global resultss
    resultss += result