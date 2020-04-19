# /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
import logging as log
import os
import io
import re
import math
import shlex
from Bio import Phylo
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import ParsimonyScorer
from heist import seqtools
from ete3 import Tree
from collections import OrderedDict
from subprocess import Popen, PIPE

"""
Hemiplasy Tool
Authors: Matt Gibson, Mark Hibbins
Indiana University
"""

def names2ints(newick):
    t = Tree(newick, format = 1)
    ntaxa = 0
    for n1 in t.iter_leaves():
        ntaxa += 1

    i = 0

    rankings = {}

    for n1 in t.iter_leaves():
        n_ancestors = len(n1.get_ancestors())
        rankings[n1.name] = n_ancestors

    s = [(k, rankings[k]) for k in sorted(rankings, key=rankings.get, reverse=False)]

    int_names = [i for i in range(1, len(s)+1)]
    rankings = {}

    for i, item in enumerate(s):
        name = item[0]
        new = int_names[i]
        newick = newick.replace(name, str(new))
        rankings[name] = new
    newick = Tree(newick, format = 1)
    newick.convert_to_ultrametric()
    return(newick.write(), rankings)



def newick2ms(newick):
    """
    Converts a Newick tree with branch lengths in coalscent units to ms-style splits
    """
    t = Tree(newick, format=1)

    ntaxa = 0
    for n1 in t.iter_leaves():
        ntaxa += 1
    i = 0
    ms_splits = []
    ms_taxa = []
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

        ms_splits.append((minDistance/2)/2)
        ms_taxa.append([max(taxa_to_collapse), min(taxa_to_collapse)])

        for node1 in t.iter_leaves():
            if node1.name == str(max(taxa_to_collapse)):
                node1.delete(preserve_branch_length=True)
        i += 1
    return(ms_splits, ms_taxa)


def splits_to_ms(splitTimes, taxa, reps, path_to_ms, y, prefix, admix=None):
    """
    Converts inputs into a call to ms
    TODO: Add introgression
    """
    nsamples = len(splitTimes) + 1
    call = (
        path_to_ms + " " + str(nsamples) + " " + str(reps) + " -T -I " + str(nsamples)
    )
    for i in range(0, nsamples):
        call += " 1"
    for x, split in enumerate(splitTimes):
        call += " -ej " + str(split) + " " + str(taxa[x][0]) + " " + str(taxa[x][1])

    if admix is not None:
        call += (
            " -es "
            + admix[0]
            + " "
            + admix[1]
            + " "
            + "1"
            + " -ej "
            + admix[0]
            + " "
            + str(nsamples + 1)
            + " "
            + admix[2]
        )

    if admix is not None:
        call += " | tail -n +4 | grep -v // > " + prefix + ".trees" + str(y) + ".tmp"
    else:
        call += " | tail -n +4 | grep -v // > " + prefix + ".trees" + str(y) + ".tmp"
    return call


def seq_gen_call(treefile, path, s, i, prefix):
    """
    Make seq-gen call.
    """
    return path + " -m HKY -l 1 -s " + str(s) + ' -wa <"' + treefile + '" > ' + prefix + '.seqs' + str(i) + '.tmp'

def print_banner():
    print(" _   _      ___ ____ _____ ")
    print("| | | | ___|_ _/ ___|_   _|")
    print("| |_| |/ _ \| |\___ \ | |  ")
    print("|  _  |  __/| | ___) || |  ")
    print("|_| |_|\___|___|____/ |_|  ")
    print("Hemiplasy Inference Simulation Tool")
    print("Version 0.3.0")
    print()
    print("Written by Mark Hibbins & Matt Gibson")
    print("Indiana University")
    print()

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
        process_ms = Popen(concatCall, stdout = PIPE, stderr=PIPE)
        #os.system(concatCall)
    else:
        log.debug("Calling ms...")
        #os.system(ms_call)
        process_ms = Popen(ms_call, shell = True)

    return(process_ms)

def call_programs_sg(ms_call, seqgencall, treefile, ntaxa):
    
        """
        Calls ms and seq-gen
        """
        log.debug("Calling seq-gen...")
        #seqgencall = shlex.split(seqgencall)
        process = Popen(seqgencall, shell = True)
        #os.system(seqgencall)
        return(process)


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
    return [c_disc_follow, c_conc_follow]


def add_branch_lengths(tree, const=1):
    """Adds dummy branch lengths to newick tree"""
    newTree = re.sub(r"(\d+)", r"\1:1", tree)
    newTree = re.sub(r"(\))", r"\1:1", newTree)
    newTree = newTree[:-2]
    return newTree


def get_min_mutations(tree, traitpattern):
    """
    Gets the minimum number of mutations required
    to explain the trait pattern without hemiplasy;
    ie. the parsimony score. Takes a newick tree 
    and trait pattern in alignment form
    """
    scorer = ParsimonyScorer()

    return scorer.get_score(tree, traitpattern)


def fitchs_alg(tree, traits):
    #tree = add_branch_lengths(tree)
    tree = Phylo.read(io.StringIO(tree), "newick")

    records = []
    for key, val in traits.items():
        if val == 0:
            records.append(SeqRecord(Seq("A", generic_dna), id=str(key)))
        elif val == 1:
            records.append(SeqRecord(Seq("T", generic_dna), id=str(key)))

    test_pattern = MultipleSeqAlignment(records)

    return get_min_mutations(tree, test_pattern)


def write_output(
    summary,
    mutation_counts_c,
    mutation_counts_d,
    reduced,
    counts,
    speciesTree,
    admix,
    traits,
    min_mutations_required,
    filename,
    reps,
    conversions,
    oldTree,
    intercept,
    coef,
    newick_internals,
    coal_internals):
    out1 = open(filename+'.txt', "w")
    out2 = open(filename+'_raw.txt', "w")

    # CALCULATE SUMMARY STATS
    #print(conversions)
    derived = []
    tree = speciesTree
    for key, val in traits.items():
        if val == 1:
            derived.append(str(key))
            tree = re.sub(r"\b%s\b" % str(key)+":", str(key) + "*:", tree)
            #tree = tree.replace(str(key)+":", (str(key) + "*:"))
    if min_mutations_required != 2:
        mix_range = list(range(2, min_mutations_required))
    else:
        mix_range = [0]
    true_hemi = 0
    mix = 0
    true_homo = 0
    for item in mutation_counts_d:
        if item[0] == 1:
            true_hemi = item[1]
        elif item[0] in mix_range:
            mix += item[1]
        elif item[0] >= min_mutations_required:
            true_homo += item[1]
    for item in mutation_counts_c:
        if item[0] >= min_mutations_required:
            true_homo += item[1]

    sum_from_introgression = counts[1]
    sum_from_species = counts[0]


    mutation_counts_comb = {}
    mutation_counts_keys = set()
    mutation_counts_cc = {}
    mutation_counts_dd = {}

    for x in mutation_counts_c:
        key = x[0]
        val = x[1]
        mutation_counts_cc[key] = val
        mutation_counts_keys.add(key)
    for x in mutation_counts_d:
        key = x[0]
        val = x[1]
        mutation_counts_dd[key] = val
        mutation_counts_keys.add(key)

    for k in mutation_counts_keys:
        if k in mutation_counts_cc.keys() and k in mutation_counts_dd.keys():
            mutation_counts_comb[k] = mutation_counts_cc[k] + mutation_counts_dd[k]
        elif k in mutation_counts_cc.keys() and k not in mutation_counts_dd.keys():
            mutation_counts_comb[k] = mutation_counts_cc[k]
        elif k not in mutation_counts_cc.keys() and k  in mutation_counts_dd.keys():
            mutation_counts_comb[k] = mutation_counts_dd[k]


    newick_internals = [str(x) for x in newick_internals]
    coal_internals = [str(x) for x in coal_internals]


    # INPUT SUMMARY
    out1.write("### INPUT SUMMARY ###\n\n")

    out1.write("Integer Code\tTaxon Name\n")
    for key, val in conversions.items():
        out1.write(str(val) + ":\t" + key + '\n')
    out1.write('\n')

    if oldTree == speciesTree:
        out1.write("The species tree (smoothed, in coalescent units) is:\n " + speciesTree + "\n\n")
    else:
        out1.write("The original species tree (smoothed, in coalescent units) is:\n " + oldTree + "\n\n")
        out1.write("The pruned species tree (smoothed, in coalescent units) is:\n " + speciesTree + "\n\n")

    out1.write("Regression intercept: " + str(intercept) + '\n')
    out1.write("Regression slope: " + str(coef) + '\n')
    out1.write("X (newick internals): " + ",".join(newick_internals) + '\n')
    out1.write("Y (coalescent internals): " + ",".join(coal_internals) + '\n')

    t = tree.replace(";", "")
    t = Phylo.read(io.StringIO(t), "newick")
    Phylo.draw_ascii(t, out1, column_width=40)

    # INPUT SUMMARY
    out1.write(
        str(len(derived))
        + " taxa have the derived state: "
        + ", ".join(derived)
        + "\n\n"
    )

    out1.write(
        "With homoplasy only, "
        + str(min_mutations_required)
        + " mutations are required to explain this trait pattern (Fitch parsimony)"
        + "\n\n"
    )

    for event in admix:
        out1.write(
            "Introgression from taxon "
            + event[1]
            + " into taxon "
            + event[2]
            + " occurs at time "
            + event[0]
            + " with probability "
            + event[3]
            + "\n"
        )
    out1.write("\n")

    out1.write(str("{:.2e}".format(reps)) + " simulations performed")

    # OUTPUT SUMMARY
    out1.write("\n\n### RESULTS ###\n\n")
    out1.write(str(sum([true_hemi, mix, true_homo])) + ' loci matched the species character states\n\n')
    out2.write(str(sum([true_hemi, mix, true_homo])) + '\n') #1#

    out1.write(
        '"True" hemiplasy (1 mutation) occurs ' + str(true_hemi) + " time(s)\n\n"
    )
    out2.write(str(true_hemi) + '\n') #2# 

    try:
        out1.write(
            "Combinations of hemiplasy and homoplasy (1 < # mutations < "
            + str(min_mutations_required)
            + ") occur "
            + str(mix)
            + " time(s)\n\n"
        )
        out2.write(str(mix) + '\n') #3#
    except:
        out1.write(
            "Combinations of hemiplasy and homoplasy (1 < # mutations < "
            + str(min_mutations_required)
            + ") occur "
            + str(0)
            + " time(s)\n\n"
        )
        out2.write(str(0) + '\n') #3#


    out1.write(
        '"True" homoplasy (>= ' + str(min_mutations_required) + ' mutations) occurs ' + str(true_homo) + " time(s)\n\n"
    )
    out2.write(str(true_homo) + '\n') #4#

    out1.write(str(summary[0]) + " loci have a discordant gene tree\n")
    out2.write(str(summary[0]) + '\n') #5#

    out1.write(
        str(summary[1] - summary[0]) + " loci are concordant with the species tree\n\n"
    )
    out2.write(str(summary[1] - summary[0]) + '\n') #6#

    out1.write(
        str(sum_from_introgression) + " loci originate from an introgressed history\n"
    )
    out2.write(str(sum_from_introgression) + '\n') #7#


    out1.write(str(sum_from_species) + " loci originate from the species history\n\n")
    out2.write(str(sum_from_species) + '\n') #8#

    # DETAILED OUTPUT
    out1.write('Distribution of mutation counts:\n\n')
    out1.write("# Mutations\t# Trees\n")
    
    out1.write("On all trees:\n")
    for key, val in mutation_counts_comb.items():
        out1.write(str(key) + '\t\t' + str(val) + '\n')
        out2.write('All' + "," + str(key) + "," + str(val) + '\n')
    out1.write("\nOn concordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in mutation_counts_c:
        out1.write(str(item[0]) + "\t\t" + str(item[1]) + "\n")
        out2.write('Conc' + "," + str(item[0]) + "," + str(item[1]) + '\n')
    out1.write("\nOn discordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in mutation_counts_d:
        out1.write(str(item[0]) + "\t\t" + str(item[1]) + "\n")
        out2.write('Disc' + "," + str(item[0]) + "," + str(item[1]) + '\n')


    if reduced is not None:
        out1.write(
            "\nOrigins of mutations leading to observed character states for hemiplasy + homoplasy cases:\n\n"
        )
        out1.write("\tTip mutation\tInternal branch mutation\tTip reversal\n")
        for key, val in reduced.items():
            val = [str(v) for v in val]
            if key in derived:
                out1.write("Taxa " + key + "\t" + "\t".join(val) + "\t" + "0" + "\n")
                out2.write("Taxa " + key + "," + ",".join(val) + "," + "0" + "\n")
            else:
                out1.write("Taxa " + key + "\t" + "\t".join(["0", val[1]]) + "\t" + val[0] + "\n")
                out2.write("Taxa " + key + "," + ",".join(["0", val[1]]) + "," + val[0] + "\n")

    out1.close()
    out2.close()


def plot_mutations(mutation_counts_c, mutation_counts_d, filename):
    """
    Plot mutation distribution with matplotlib
    """

    mutation_counts_comb = {}
    mutation_counts_keys = set()
    mutation_counts_cc = {}
    mutation_counts_dd = {}

    for x in mutation_counts_c:
        key = x[0]
        val = x[1]
        mutation_counts_cc[key] = val
        mutation_counts_keys.add(key)
    for x in mutation_counts_d:
        key = x[0]
        val = x[1]
        mutation_counts_dd[key] = val
        mutation_counts_keys.add(key)

    for k in mutation_counts_keys:
        if k in mutation_counts_cc.keys() and k in mutation_counts_dd.keys():
            mutation_counts_comb[k] = mutation_counts_cc[k] + mutation_counts_dd[k]
        elif k in mutation_counts_cc.keys() and k not in mutation_counts_dd.keys():
            mutation_counts_comb[k] = mutation_counts_cc[k]
        elif k not in mutation_counts_cc.keys() and k  in mutation_counts_dd.keys():
            mutation_counts_comb[k] = mutation_counts_dd[k]



    objs_comb = mutation_counts_comb.keys()
    conc_dic = {}
    disc_dic = {}

    objs = set(objs_comb)
    x = np.array(list(range(1, max(objs) + 1)))
    y1 = []
    width = 0.2


    for i in x:
        if i in mutation_counts_comb.keys():
            y1.append(mutation_counts_comb[i])
        else:
            y1.append(0)

    labels = []
    for o in x:
        labels.append(str(o))
    _, ax = plt.subplots()
    p1 = ax.bar(x, y1, width, color="#484041")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    plt.ylabel("Count")
    plt.xlabel("# Mutations")
    plt.savefig(filename + ".dist.png", dpi=50)


def subs2coal(newick_string):
        
        '''
        Takes a newick string with the nodes labelled with concordance factors,
        and returns the same string with branch lengths converted to coalescent
        units.
        '''

        scfs = re.findall("\)(.*?)\:", newick_string) #regex to get concordance factors
        scfs = list(filter(None, scfs)) #removes empty values (root has no scf)
        coal_internals = []

        for i in range(len(scfs)):

                if (float(scfs[i])/100) < 1 and float(scfs[i]) != 1: #catch for missing data / values of 1
                        coal_estimate = -1*(math.log(3/2) + math.log(1 - float(scfs[i])/100))
                        if coal_estimate > 0:
                                coal_internals.append(coal_estimate)
                        else:
                                coal_internals.append(0.01)
                else:
                        coal_internals.append(np.NaN)

        coal_internals = [float(i) for i in coal_internals]
        newick_internals = []
        internal_sections = []
        sections = newick_string.split(':')

        for i in range(len(sections)):
                if len(sections[i].split(")")) > 1 and sections[i].split(")")[1] in scfs:
                                newick_internals.append(newick_string.split(':')[i+1].split(',')[0].split(")")[0])
                                internal_sections.append(sections[i+1])

        
        newick_internals = [float(i) for i in newick_internals]
         
        def branch_regression(newick_branches, coal_internals):

                '''
                Returns the regression coefficient of the coalescent branch lengths
                on the branch lengths in the newick string
                '''

                for i in range(len(newick_branches)): #drops missing data
                        if np.isnan(newick_branches[i]) or np.isnan(coal_internals[i]):
                                del newick_branches[i]
                                del coal_internals[i]

                intercept, slope = poly.polyfit(newick_branches, coal_internals, 1)

                return intercept, slope

        intercept, coef = branch_regression(newick_internals, coal_internals)
        n = newick_internals
        c = coal_internals
        tip_sections = []

        for i in range(len(sections)): #gets the sections with tip lengths
                if sections[i] not in internal_sections:
                        tip_sections.append(sections[i])

        newick_tips = []

        for i in range(len(tip_sections)): #gets the tip lengths
                if len(tip_sections[i].split(',')) > 1:
                        newick_tips.append(tip_sections[i].split(',')[0])
                elif len(tip_sections[i].split(')')) > 1:
                        newick_tips.append(tip_sections[i].split(')')[0])
        
        newick_tips = [float(i) for i in newick_tips]

        coal_tips = [(newick_tips[i]*coef + intercept) for i in range(len(newick_tips))]

        for i in range(len(coal_tips)):
                if coal_tips[i] <= 0:
                        coal_tips[i] = 0.01
                else:
                        coal_tips[i] = coal_tips[i]

        coal_internals = [(newick_internals[i]*coef + intercept) if np.isnan(coal_internals[i]) else coal_internals[i] for i in range(len(newick_internals))]

        lengths = re.findall("\d+\.\d+", newick_string)

        lengths = [lengths[i] for i in range(len(lengths)) if float(lengths[i]) < 1]

        coal_lengths = []

        for i in range(len(lengths)):
                if newick_internals.count(lengths[i]) > 1 or newick_tips.count(lengths[i]) > 1:
                        sys.exit('Error: Duplicate branch lengths')
                elif float(lengths[i]) in newick_internals:
                        coal_lengths.append(coal_internals[newick_internals.index(float(lengths[i]))])
                elif float(lengths[i]) in newick_tips:
                        coal_lengths.append(coal_tips[newick_tips.index(float(lengths[i]))])
                elif float(lengths[i]) == 0: #deals with roots of length 0 in smoothed trees
                        coal_lengths.append(float(0))
        
        coal_lengths = [str(coal_lengths[i]) for i in range(len(coal_lengths))]
       
        coal_newick_string = newick_string

        for i in range(len(lengths)):
             coal_newick_string = coal_newick_string.replace(str(lengths[i]), str(coal_lengths[i]))
    
        scfs = [float(x) for x in scfs]
        scfs2 = []
        for val in scfs:
            if val.is_integer():
                scfs2.append(str(val) + ".0")
            else:
                scfs2.append(str(val))


        for i in range(len(scfs)):
                coal_newick_string = coal_newick_string.replace(str(scfs[i]), '')    
        return(coal_newick_string, Tree(coal_newick_string, format=1), intercept, coef, n, c)


def readInput(file):
    f = open(file, "r")
    tree = ""
    derived = []
    admix = []
    outgroup=None
    cnt = 0
    treeType = 'ml'

    for i, line in enumerate(f):
        if line.startswith("begin trees"):
            cnt = 1
            continue
        if line.startswith("begin hemiplasytool"):
            cnt = 2
            continue
        if (len(line) <= 1) or (line.startswith("end")) or (line.startswith("#")):
            continue


        elif cnt == 1:
            tree = line.replace("\n", "").split('=')[1][1:]

        elif cnt == 2:
            #simulation parameters
            l = line.replace("\n", "").split('=')
            if l[0] == "set derived taxon":
                derived.append(l[1])
            elif l[0] == 'set outgroup taxon':
                outgroup = l[1]
            elif l[0].startswith('set introgression'):
                l = line.replace('\n','').split(' ')
                sp1 = l[2].split('=')[1]
                sp2 = l[3].split('=')[1]
                strength = l[4].split('=')[1]
                time = l[5].split('=')[1]
                admix.append([time,sp1,sp2,strength])
            elif l[0].startswith('set type coal'):
                treeType = 'coal'
    
            
    return(tree, derived, admix, outgroup, treeType)


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

    # print('\nDerived mutation inheritance patterns for trees with fewer\
    #     mutations than derived taxa:\n')
    # print('\tTerm\tInherited from anc node')
    # for key, val in reduced.items():
    #    val = [str(v) for v in val]
    #    print('Taxa ' + key + '\t' + '\t'.join(val))
    return reduced


def update_count(tree, dic):
    for key, val in dic.items():
        if seqtools.compareToSpecies(key, tree):
            dic[key] += 1
    return dic


def write_unique_trees(focal_trees, filename, traits):
    unique = []
    counts = {}
    out1 = open(filename+'.txt', "a")
    outTrees = open(filename+'.trees', 'w')
    for i, tree in enumerate(focal_trees):
        outTrees.write(tree + '\n')
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
    out1.write("\n### OBSERVED GENE TREES ###\n\n")
    for tree in unique:
        for key, val in traits.items():
            if val == 1:
                tree = re.sub(r"\b%s\b" % str(key)+":", str(key) + "*:", tree)
        t = tree
        t = t.replace(";", "")
        t = Phylo.read(io.StringIO(t), "newick")
        Phylo.draw_ascii(t, out1, column_width=40)
        for key, val in counts.items():

            if seqtools.compareToSpecies(key, tree):
                out1.write("This topology occured " + str(val) + " time(s)\n")

    out1.close()
    outTrees.close()


def prune_tree(tree, derived, outgroup):
    t = Tree(tree, format=1)
    for leaf in t:
        if leaf.name in derived:
            leaf.add_features(derived='1')
        else:
            leaf.add_features(derived='0')
    ns = []
    for node in t.get_monophyletic(values=["1"], target_attr="derived"):
       ns.append(node)
    sub = t.get_common_ancestor(ns[0], ns[1], ns[2], ns[3])
    tokeep = list(sub.get_leaf_names())
    tokeep.append(outgroup)
    t.prune(tokeep)
    return(t.write(), t)
    
