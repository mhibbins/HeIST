from ete3 import Tree
from itertools import zip_longest

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return(zip_longest(*args, fillvalue=fillvalue))

def cluster(d):
    clusters = {}
    for key, val in d.items():
        clusters.setdefault(val, []).append(key)
    return(clusters)

def checkEqual(lst):
   return(lst[1:] == lst[:-1])


def readSeqs(seqs, ntaxa, speciesPattern=None):
    indices = []
    c = cluster(speciesPattern)
    shouldMatch1 = c['0']
    shouldMatch2 = c['1']

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
    focal_trees = []
    all_trees = []
    trees = open(treefile, 'r')
    for i, line in enumerate(trees):
        l = line.replace('\n','')
        if len(l) > 3:
            all_trees.append(l)
    trees.close()
    for i, tree in enumerate(all_trees):
        if i in matchlist:
            focal_trees.append(tree)
    return(focal_trees)
            
def compareToSpecies(speciesTree, geneTree):
    speciesTree = Tree(speciesTree)
    geneTree = Tree(geneTree)
    r = speciesTree.compare(geneTree)['rf']
    if r == 0.0:
        return(True)
    else:
        print(geneTree)
        return(False)

def propDiscordant(focal_trees, species_tree):
    countDis = 0

    for tree in focal_trees:
        if not compareToSpecies(species_tree, tree):
            countDis += 1
    return(countDis, len(focal_trees), countDis/len(focal_trees))

