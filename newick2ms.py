#Tester script for converting newick trees to ms inputs

from ete3 import Tree

newick = "(1:7.0,(2:4.0,((6:1.5,5:1.5):1.5,(4:1.0,3:1.0):2.0):1.0):3.0);"
        
newick2 = "(1:10.0,(2:8.0,((3:2,(4:1,5:1):1):5.0,6:7.0):1.0):2.0);"
t = Tree(newick2, format=1)
    
ntaxa = 0
for n1 in t.iter_leaves():
    ntaxa += 1


i = 0

ms_string = ""

while i < ntaxa+1:
    distances = []
    taxa = []

    for n1 in t.iter_leaves():
        for n2 in t.iter_leaves():
            n1n = n1.name
            n2n = n2.name
            d = n1.get_distance(n2)
            if n1n != n2n:
                distances.append(d)
                taxa.append(n1n + ' ' + n2n)

    if len(distances) == 0:
        break
    minDistance_idx = distances.index(min(distances))
    minDistance = min(distances)

    taxa_to_collapse = taxa[minDistance_idx].split(' ')
    taxa_to_collapse = [int(x) for x in taxa_to_collapse]
    ms_string = ms_string + str(max(taxa_to_collapse)) + ' ' + str(min(taxa_to_collapse)) + ' ' + str(minDistance/2) + '\n'

    for node1 in t.iter_leaves():
        if node1.name == str(max(taxa_to_collapse)):
            node1.delete(preserve_branch_length=True)
    i += 1
print(ms_string)