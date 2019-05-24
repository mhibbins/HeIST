# /usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import logging as log
import os
import io
import re
from Bio import Phylo
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import ParsimonyScorer
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
        call += " | tail -n +4 | grep -v // > trees" + str(r) + ".tmp"
    else:
        call += " | tail -n +4 | grep -v // > trees.tmp"
    return call


def seq_gen_call(treefile, path, s=0.05):
    """
    Make seq-gen call.
    """
    return path + " -m HKY -l 1 -s " + str(s) + ' -wa <"' + treefile + '" > seqs.tmp'


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
    tree = add_branch_lengths(tree)
    tree = Phylo.read(io.StringIO(tree), "newick")

    records = []
    for key, val in traits.items():
        if val == "0":
            records.append(SeqRecord(Seq("A", generic_dna), id=key))
        elif val == "1":
            records.append(SeqRecord(Seq("T", generic_dna), id=key))

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
):
    out1 = open(filename, "w")

    # CALCULATE SUMMARY STATS
    derived = []
    tree = speciesTree
    for key, val in traits.items():
        if val == "1":
            derived.append(key)
            tree = tree.replace(key, (key + "*"))
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

    sum_from_introgression = 0
    sum_from_species = 0
    for i, topology_count in enumerate(counts):
        if i != len(counts) - 1:
            sum_from_introgression += topology_count
        else:
            sum_from_species = topology_count

    # INPUT SUMMARY
    out1.write("### INPUT SUMMARY ###\n\n")
    out1.write("The species tree is " + speciesTree + "\n\n")
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
        "The minimum number of mutations required to explain this trait pattern is "
        + str(min_mutations_required)
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

    # OUTPUT SUMMARY
    out1.write("\n### OUTPUT SUMMARY ###\n\n")
    out1.write(
        '"True" hemiplasy (1 mutation) occurs ' + str(true_hemi) + " time(s)\n\n"
    )
    if mix_range != [0]:
        out1.write(
            "Combinations of hemiplasy and homoplasy (1 < # mutations < "
            + str(min_mutations_required)
            + ") occur "
            + str(mix)
            + " time(s)\n\n"
        )
    out1.write(
        '"True" homoplasy (>= 3 mutations) occurs ' + str(true_homo) + " time(s)\n\n"
    )
    out1.write(str(summary[0]) + " loci have a discordant gene tree\n")
    out1.write(
        str(summary[1] - summary[0]) + " loci are concordant with the species tree\n\n"
    )
    out1.write(
        str(sum_from_introgression) + " loci originate from an introgressed history\n"
    )
    out1.write(str(sum_from_species) + " loci originate from the species history\n\n")
    if reduced is not None:
        out1.write("In cases with combinations of hemiplasy and homoplasy:\n\n")
        for key, val in reduced.items():
            val = [str(v) for v in val]
            out1.write(
                "Taxon "
                + key
                + " mutated to the derived state "
                + val[0]
                + " time(s), and inherited it from an ancestral population "
                + val[1]
                + " time(s)\n"
            )

    # DETAILED OUTPUT
    out1.write("\n\n### DETAILED OUTPUT ###\n\n")
    out1.write("On concordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in mutation_counts_c:
        out1.write(str(item[0]) + "\t\t" + str(item[1]) + "\n")
    out1.write("\nOn discordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in mutation_counts_d:
        out1.write(str(item[0]) + "\t\t" + str(item[1]) + "\n")

    if reduced is not None:
        out1.write(
            "\nDerived mutation inheritance patterns for trees with fewer mutations than derived taxa:\n"
        )
        out1.write("\tTerm\tInherited from anc node\n")
        for key, val in reduced.items():
            val = [str(v) for v in val]
            out1.write("Taxa " + key + "\t" + "\t".join(val) + "\n")

    out1.write("\nOf the replicates that follow species site pattern:\n")
    out1.write(
        str(summary[0])
        + " were discordant\n"
        + str(summary[1] - summary[0])
        + " were concordant\n\n"
    )

    for i, topology_count in enumerate(counts):
        if i != len(counts) - 1:
            out1.write(
                str(topology_count)
                + " replicates matching the species pattern were from introgression tree "
                + str(i + 1)
                + "\n"
            )
        else:
            out1.write(
                str(topology_count)
                + " replicates matching the species pattern were from the species tree\n"
            )

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
    x = np.array(list(range(1, max(objs) + 1)))
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
    p1 = ax.bar(x, y1, width, color="#484041")
    p2 = ax.bar(x + width, y2, width, color="#70ee9c")
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(labels)
    plt.ylabel("Count")
    plt.xlabel("# Mutations")
    ax.legend((p1[0], p2[0]), ("Concordant trees", "Discordant trees"))
    plt.savefig("mutation_dist.png", dpi=250)


def readInput(file):
    f = open(file, "r")
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

    return (splits, taxa, traits, tree, admix)


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


def write_unique_trees(focal_trees, filename):
    unique = []
    counts = {}
    out1 = open(filename, "a")

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
    out1.write("\n### OBSERVED GENE TREES ###\n\n")
    for tree in unique:
        t = tree.replace(";", "")
        t = Phylo.read(io.StringIO(tree), "newick")
        Phylo.draw_ascii(t, out1, column_width=40)
        for key, val in counts.items():

            if seqtools.compareToSpecies(key, tree):
                out1.write("This topology occured " + str(val) + " time(s)\n")

    out1.close()
