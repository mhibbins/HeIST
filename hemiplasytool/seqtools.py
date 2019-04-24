"""
Scripts that do the parsing of gene trees and sequences.
"""

from ete3 import Tree
from itertools import zip_longest
from Bio.Phylo.Consensus import _BitString
from Bio import Phylo
import io

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
            print(pattern)
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

def compareToSpecies(tree1, tree2):
    tree1 = tree1.replace(";", '')
    tree2 = tree2.replace(';', '')
    tree1 = Phylo.read(io.StringIO(tree1), "newick")
    tree2 = Phylo.read(io.StringIO(tree2), "newick")
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False

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
    Determines the proportion of focal_trees (which have the same site pattern as the
    species tree) which are discordant (i.e. have a different topology)
    """
    countDis = 0

    for tree in focal_trees:
        if not compareToSpecies(species_tree, tree):
            countDis += 1
    try:
        return(countDis, len(focal_trees), countDis/len(focal_trees))
    except:
        return(countDis, len(focal_trees), 0.0)


