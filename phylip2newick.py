import re
import numpy as np
file1 = "./test_phylip.txt"
ntaxa = 15

def replace(inputcase, case):
    print(inputcase)
    print(case)
    new = [case[0],case[1]]
    inputcase.remove(case[0])
    for i,item in enumerate(inputcase):
        if item == case[1]:
            inputcase[i] = new
    print(inputcase)
    return(inputcase)

def recursive_len(item):
    if type(item) == list:
        return sum(recursive_len(subitem) for subitem in item)
    else:
        return 1

def parse_seqgen(seqfile, ntaxa):
    '''
    Parses seq-gen file. Returns a list of ordered 
    taxon-allele pairs for each tree.   
    '''
    lines = []

    with open(seqfile) as seqs:
        for line in seqs:
            if re.match(r'\w', line): #if line starts with character
                l = line.split()
                lines.append(l[0])

    #lines = [lines[i].replace('\t', ' ') for i in range(len(lines))] #replaces tabs with space
	
    trees = [lines[i:i+(ntaxa*2)-1] for i in range(0, len(lines), (ntaxa*2)-1)] #splits into gene trees
    return trees

tree = parse_seqgen(file1,ntaxa)[0]

print(tree)

possible = [str(i) for i in list(range(1,ntaxa+1))]


complete = False
done = False
x = 0

while(complete == False):
    #Merge sister taxa into clades
    while(done == False):
        for i,taxa in enumerate(tree):
            x += 1
            if x > len(tree)+50:
                done = True
            if i != 0:
                prev = tree[i-1]
                case = (taxa,prev)
                #print(case)
                if (case[0] in possible) and (case[1] in possible):
                    tree = replace(tree, case)
                    break
        print(tree)

    for i,item in enumerate(tree):
        if not isinstance(item, list):
            tree[i] = [tree[i]]
    print(tree)
    


    done = False
    x = 0

    #Merge pairs of sisters (seperated by one ancestral node)
    while(done == False):
        for i,taxa in enumerate(tree):
            x += 1
            if x > len(tree):
                done = True
            if (i != 0) and (i < len(tree)):
                try:
                    one = tree[i-1]
                    two = tree[i]
                    three = tree[i+1]
                except IndexError:
                    done = True
                    break
                if (recursive_len(one) > 1) and (recursive_len(three) > 1) and (two[0] not in possible):
                    tree.remove(two)
                    tree.remove(three)
                    tree[i-1] = [one] + [three]
                    break
                elif (one[0] in possible) and (recursive_len(three) > 1) and (len(two) == 1) and (two[0] not in possible):
                    tree.remove(two)
                    tree.remove(three)
                    tree[i-1] = [[one[0]] + [three]]
                    break
    print(tree)



    done = False
    x = 0
    #Merge clades seperated by two ancestral nodes
    while(done == False):
        for i,taxa in enumerate(tree):
            x += 1
            if x > len(tree):
                done = True
            if (i != 0) and (i < len(tree)):
                try:
                    one = tree[i-1]
                    two = tree[i]
                    three = tree[i+1]
                    four = tree[i+2]
                except IndexError:
                    done = True
                    break

                #
                if (one[0] in possible) and (two[0] not in possible) and (three[0] not in possible) and (len(four) > 1):
                    tree.remove(two)
                    tree.remove(three)
                    tree.remove(four)
                    tree[i-1] = [[one[0]]+[four]]
                    break
                elif (len(one) > 1) and (two[0] not in possible) and (three[0] not in possible) and (four[0] in possible):
                    tree.remove(two)
                    tree.remove(three)
                    tree.remove(one)
                    tree[i-1] = [[one] + [four[0]]]

    print(tree)
    if (recursive_len(tree) == ntaxa):
        complete = True
        break

    if tree[0][0] not in possible:
        tree.remove(tree[0])
        if (recursive_len(tree) == ntaxa):
            complete = True
            break
    else:
        print("NO OUTGROUP")
        exit(0)
    


    print(tree)
print(tree)