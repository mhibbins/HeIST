# /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import logging as log
import os
import io
from Bio import Phylo
from hemiplasytool import seqtools


"""
Hemiplasy Tool
Authors: Matt Gibson, Mark Hibbins
Indiana University
"""


def splits_to_ms(splitTimes, taxa, reps, path_to_ms, admix=None, r=None):
    """
    Converts inputs into a call to ms
    TODO: Add introgression
    """
    nsamples = len(splitTimes)+1
    call = path_to_ms + ' ' + str(nsamples) + ' ' + \
        str(reps) + ' -T -I ' + str(nsamples)
    for i in range(0, nsamples):
        call += ' 1'
    for x, split in enumerate(splitTimes):
        call += ' -ej ' + str(split) + ' ' + str(taxa[x][0]) + ' ' +\
            str(taxa[x][1])

    if admix is not None:
        call += ' -es ' + admix[0] + ' ' + admix[1] + ' ' + '1' + ' -ej ' +\
            admix[0] + ' ' + str(nsamples+1) + ' ' + admix[2]

    if admix is not None:
        call += " | tail -n +4 | grep -v // > trees" + str(r) + ".tmp"
    else:
        call += " | tail -n +4 | grep -v // > trees.tmp"
    return(call)


def seq_gen_call(treefile, path, s=0.05):
    """
    Make seq-gen call.
    """
    return(path + ' -m HKY -l 1 -s ' + str(s) + ' -wa <"' +
                                                treefile + '" > seqs.tmp')


def call_programs(ms_call, seqgencall, treefile, ntaxa):
    """
    Calls ms and seq-gen
    """
    if type(ms_call) is list:
        for call in ms_call:
            log.debug("Calling ms...")
            os.system(call)

        concatCall = "cat "
        for i, call in enumerate(ms_call):
            concatCall += "trees" + str(i) + ".tmp "
        concatCall += "> trees.tmp"
        os.system(concatCall)
    else:
        log.debug("Calling ms...")
        os.system(ms_call)

    log.debug("Calling seq-gen...")
    os.system(seqgencall)


def cleanup():
    """Remove gene trees and sequences files. For use between batches."""
    os.system("rm trees.tmp")
    os.system("rm seqs.tmp")
    os.system("rm focaltrees.tmp")


def summarize(results):
    """Summarizes simulations from multiple batches"""
    c_disc_follow = 0
    c_conc_follow = 0
    for key, val in results.items():
        c_disc_follow += val[0]
        c_conc_follow += val[1]
    return([c_disc_follow, c_conc_follow])


def write_output(summary, mutation_counts_c, mutation_counts_d, reduced, filename):
    out1 = open(filename, 'w')

    out1.write("On concordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in mutation_counts_c:
        out1.write(str(item[0]) + '\t\t' + str(item[1]) + '\n')
    out1.write("\nOn discordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in mutation_counts_d:
        out1.write(str(item[0]) + '\t\t' + str(item[1]) + '\n')

    if reduced is not None:
        out1.write('\nDerived mutation inheritance patterns for trees with fewer mutations than derived taxa:\n')
        out1.write('\tTerm\tInherited from anc node\n')
        for key, val in reduced.items():
            val = [str(v) for v in val]
            out1.write('Taxa ' + key + '\t' + '\t'.join(val) + '\n')
    
    out1.write("\nOf the replicates that follow species site pattern:\n")
    out1.write(str(summary[0]) + " were discordant\n" + 
    str(summary[1]-summary[0]) + " were concordant\n")


    out1.close()


def plot_mutations(results_c, results_d):
    """
    Plot mutation distribution with matplotlib
    """
    objs_c = [i[0] for i in results_c]
    objs_d = [i[0] for i in results_d]
    conc_dic = {}
    disc_dic = {}
    objs = objs_c + objs_d
    objs = set(objs)
    x = np.array(list(range(1, max(objs)+1)))
    y1 = []
    y2 = []
    width = 0.2
    for v in results_c:
        conc_dic[v[0]] = v[1]
    for v in results_d:
        disc_dic[v[0]] = v[1]

    for i in x:
        if i in conc_dic.keys():
            y1.append(conc_dic[i])
        else:
            y1.append(0)
        if i in disc_dic.keys():
            y2.append(disc_dic[i])
        else:
            y2.append(0)

    labels = []
    for o in x:
        labels.append(str(o))
    _, ax = plt.subplots()
    p1 = ax.bar(x, y1, width, color='#484041')
    p2 = ax.bar(x+width, y2, width, color='#70ee9c')
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(labels)
    plt.ylabel('Count')
    plt.xlabel('# Mutations')
    ax.legend((p1[0], p2[0]), ('Concordant trees', 'Discordant trees'))
    plt.savefig('mutation_dist.png', dpi=250)


def readInput(file):
    f = open(file, 'r')
    splits = []
    taxa = []
    traits = {}
    admix = []

    cnt = 0
    for i, line in enumerate(f):
        if line.startswith("#"):
            cnt += 1
            continue
        if len(line) <= 1:
            continue

        elif cnt == 1:
            l = line.replace("\n", "").split()
            splits.append(float(l[0]))
            taxa.append((int(l[1]), int(l[2])))

        elif cnt == 2:
            l = line.replace("\n", "").split()
            traits[l[0]] = l[1]

        elif cnt == 3:
            tree = line.replace("\n", "")
        elif cnt == 4:
            l = line.replace("\n", "").split()
            admix.append(l)

    return(splits, taxa, traits, tree, admix)


def summarize_inherited(inherited):
    reduced = {}
    for event in inherited:
        if event[0] not in reduced.keys():
            if event[1] == 1:
                reduced[event[0]] = [1, 0]
            elif event[1] == 0:
                reduced[event[0]] = [0, 1]
        elif event[0] in reduced.keys():
            if event[1] == 1:
                reduced[event[0]][0] += 1
            elif event[1] == 0:
                reduced[event[0]][1] += 1

    print('\nDerived mutation inheritance patterns for trees with fewer\
         mutations than derived taxa:\n')
    print('\tTerm\tInherited from anc node')
    for key, val in reduced.items():
        val = [str(v) for v in val]
        print('Taxa ' + key + '\t' + '\t'.join(val))
    return(reduced)


def update_count(tree, dic):
    for key, val in dic.items():
        if seqtools.compareToSpecies(key, tree):
            dic[key] += 1
    return(dic)


def write_unique_trees(focal_trees, filename):
    unique = []
    counts = {}
    out1 = open(filename, 'a')

    for i, tree in enumerate(focal_trees):
        if i == 0:
            unique.append(tree)
            counts[tree] = 1
        else:
            uniq = True
            for tree2 in unique:
                if seqtools.compareToSpecies(tree, tree2):
                    uniq = False
                    counts = update_count(tree, counts)
            if uniq is True:
                unique.append(tree)
                counts[tree] = 1

    for tree in unique:
        t = tree.replace(";", '')
        t = Phylo.read(io.StringIO(tree), "newick")
        Phylo.draw_ascii(t, out1, column_width=40)
        for key, val in counts.items():

            if seqtools.compareToSpecies(key, tree):
                out1.write("This topology occured " + str(val) + " time(s)\n")

    out1.close()