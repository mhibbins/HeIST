# /usr/bin/python3

"""
HemiplasyTool
Authors: Matt Gibson, Mark Hibbins
Indiana University
"""

import argparse
import time
import sys
import os
import logging as log
import subprocess
from heist import hemiplasytool
from heist import seqtools
from Bio import Phylo
from ete3 import Tree

def newick2ms(*args):
    parser = argparse.ArgumentParser(
        description="Tool for converting \
                        a newick string to ms-style splits. \
                            Note that this only makes sense if the input\
                                tree is in coalescent units."
    )

    parser.add_argument(
        "input",
        metavar="input",
        help="Input newick string file"
    )
    args = parser.parse_args()

    newick = open(args.input, 'r').read()
    
    tree, conversions = hemiplasytool.names2ints(newick)

    splits, taxa = hemiplasytool.newick2ms(tree)

    call = ""
    for x, split in enumerate(splits):
        call += " -ej " + str(split) + " " + str(taxa[x][0]) + " " + str(taxa[x][1])
    print("Split times flags for ms: ")
    print(call)
    print()

    print("Time\tTaxa1\tTaxa2")
    for i, split in enumerate(splits):
        print(str(split) + '\t' + str(taxa[i][0]) + '\t' + str(taxa[i][1]))

    print('\n')

    print("Original code\tms code")
    for key, val in conversions.items():
        print(str(key) + '\t' + str(val))


def subs2coal(*args):
    parser = argparse.ArgumentParser(
        description="Tool for converting \
                        a newick string with branch lengths in subs/site\
                            to a neewick string with branch lengths in\
                                 coalescent units. Input requires gene \
                                     or site-concordancee factors as \
                                         branch labels"
    )

    parser.add_argument(
        "input",
        metavar="input",
        help="Input newick string file"
    )
    args = parser.parse_args()

    newick = open(args.input, 'r').read()

    treeSp,t = hemiplasytool.subs2coal(newick)

def heistMerge(*args):
    parser = argparse.ArgumentParser(
        description="Merge output files from multiple HeiST runs. \
            Useful for simulating large trees by running multiple batch jobs.")

    parser.add_argument(
        "-d", help="Merge all files in a directory", action='store_true'
    )
    parser.add_argument("inputs", nargs="*", help = "Prefixes of output files to merge or a directory (supply -d flag as well)")
    args = parser.parse_args()
    files = args.inputs
    
    if args.d == True:
        direc = str(files[0])
        result = subprocess.run('echo ' + str(direc) + '/*_raw.txt', stdout=subprocess.PIPE, shell=True)
        
        files2 = result.stdout.split()
        files3 = []
        for file in files2:
            file = str(file.decode("utf-8"))
            files3.append(file.replace('_raw.txt', ''))
    else:
        files3 = files
    #Get info that is same across all runs

    head = ""
    file1 = open(files3[0] + ".txt")
    for line in file1:
        if not line.startswith('### RESULTS ###'):
            head += line
        if line.startswith('### RESULTS ###'):
            break
    
    data = {"#1": 0, "#2": 0, "#3": 0, "#4": 0, "#5": 0,
            "#6": 0, "#7": 0, "#8": 0}
    allT = {}
    discT = {}
    concT = {}
    taxaT_1 = {}
    taxaT_2 = {}
    taxaT_3 = {}

    for file in files3:
        f = open(file + "_raw.txt")
        for i, line in enumerate(f):
            if i == 0:
                data["#1"] += int(line)
            elif i == 1:
                data["#2"] += int(line)
            elif i == 2:
                data["#3"] += int(line)
            elif i == 3:
                data["#4"] += int(line)
            elif i == 4:
                data["#5"] += int(line)
            elif i == 5:
                data["#6"] += int(line)
            elif i == 6:
                data["#7"] += int(line)
            elif i == 7:
                data["#8"] += int(line)
            else:
                l = line.replace('\n', '').split(',')
                if l[0] == 'All':
                    if l[1] in allT.keys():
                        allT[l[1]] += int(l[2])
                    else:
                        allT[l[1]] = int(l[2])
                elif l[0] == 'Disc':
                    if l[1] in discT.keys():
                        discT[l[1]] += int(l[2])
                    else:
                        discT[l[1]] = int(l[2])
                elif l[0] == 'Conc':
                    if l[1] in concT.keys():
                        concT[l[1]] += int(l[2])
                    else:
                        concT[l[1]] = int(l[2])
                elif l[0].startswith("Taxa"):
                    if l[0] in taxaT_1.keys():
                        taxaT_1[l[0]] += int(l[1])
                        taxaT_2[l[0]] += int(l[2])
                        taxaT_3[l[0]] += int(l[3])
                    else:
                        taxaT_1[l[0]] = int(l[1])
                        taxaT_2[l[0]] = int(l[2])
                        taxaT_3[l[0]] = int(l[3])
    print(head)

    print("### RESULTS ###")
    print(str(data["#1"]) + ' loci matched the species character states\n')
    print('"True" hemiplasy (1 mutation) occurs ' + str(data["#2"]) + " time(s)\n")
    print("Combinations of hemiplasy and homoplasy (1 < # mutations < "
            + "____"
            + ") occur "
            + str(data["#3"])
            + " time(s)\n"
        )
    print(
        '"True" homoplasy (>= ' + "______" + ' mutations) occurs ' + str(data["#4"]) + " time(s)\n"
    )
    
    print(str(data["#5"]) + " loci have a discordant gene tree\n")
    print(str(data["#6"]) + " loci are concordant with the species tree\n")
    print(str(data["#7"]) + " loci originate from an introgressed history\n")
    print(str(data["#8"]) + " loci originate from the species history\n")


    print('Distribution of mutation counts:\n')
    print("# Mutations\t# Trees")

    print("On all trees:")
    for key, val in allT.items():
        print(str(key) + '\t\t' + str(val))
    print("\nOn concordant trees:")
    for key, val in concT.items():
        print(str(key) + '\t\t' + str(val))
    print("\nOn discordant trees:")
    for key, val in discT.items():
        print(str(key) + '\t\t' + str(val))
    print()

    print("\nOrigins of mutations leading to observed character states for hemiplasy + homoplasy cases:\n")
    print("\tTip mutation\tInternal branch mutation\tTip reversal")
    for key, val in taxaT_1.items():
        print(key + "\t" + str(val) + "\t" + str(taxaT_2[key]) + "\t" + str(taxaT_3[key]))

    tree_files = [x + ".trees" for x in files3]
    catcall = "cat "
    for t in tree_files:
        catcall += t + " "
    catcall += "> merged_trees.trees"
    os.system(catcall)

def main(*args):
    start = time.time()
    hemiplasytool.print_banner()

    parser = argparse.ArgumentParser(
        description="Tool for characterising \
                        hemiplasy given traits mapped onto a species tree"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Enable debugging messages to be displayed",
        action="store_true",
    )
    parser.add_argument(
        "input",
        metavar="input",
        help="Input NEXUS file",
    )
    parser.add_argument(
        "-n",
        "--replicates",
        metavar="",
        help="Number of replicates per batch",
        default=1000000,
    )
    parser.add_argument(
        "-t", "--threads", metavar="", help="Number of threads for simulations", default=16
    )
    parser.add_argument(
        "-p", "--mspath", metavar="", help="Path to ms (if not in user path)", default="ms"
    )
    parser.add_argument(
        "-g", "--seqgenpath", metavar="", help="Path to seq-gen (if not in user path)", default="seq-gen"
    )
    parser.add_argument(
        "-s",
        "--mutationrate",
        metavar="",
        help="Seq-gen mutation rate (default 0.05)",
        default=0.05,
    )
    parser.add_argument("-o", "--outputdir", metavar="", help="Output directory/prefix")

    args = parser.parse_args()

    # Setup ###################
    log.basicConfig(level=log.DEBUG)
    logger = log.getLogger()
    if args.verbose:
        logger.disabled = False
    else:
        logger.disabled = True
    mpl_logger = log.getLogger("matplotlib")
    mpl_logger.setLevel(log.WARNING)
    ##########################
    

    # Read input file
    log.debug("Reading input file...")
    treeSp, derived, admix, outgroup, type = hemiplasytool.readInput(args.input)
    tmp1 = Tree(treeSp, format = 1)
    tmp1.convert_to_ultrametric()

    if type != 'coal':
        # Convert ML tree to a coalescent tree based on GCFs
        treeSp,t,intercept,coef,newick_internals,coal_internals = hemiplasytool.subs2coal(treeSp)
        original_tree = [treeSp, t]
    # Tree pruning
    if outgroup != None:
        log.debug("Pruning tree...")
        # Prune tree
        treeSp,t = hemiplasytool.prune_tree(treeSp, derived, outgroup)


    taxalist = [i.name for i in t.iter_leaves()]
    [i.name for i in t.iter_leaves()]

    # Convert coalescent tree to ms splits
    treeSp, conversions = hemiplasytool.names2ints(treeSp)
    original_tree[0], tmp = hemiplasytool.names2ints(original_tree[0])
    
    # Convert newick tree to ms splits
    splits, taxa = hemiplasytool.newick2ms(treeSp)
    traits = {}
    for i in taxalist:
        if i in derived:
            traits[conversions[i]] = 1
        else:
            traits[conversions[i]] = 0
    
    
    # Make program calls
    threads = int(args.threads)
    reps = int(args.replicates)

    breaks = [] #WHAT TO DO WITH THIS?????





    # Begin batches
    taxalist = []
    for s in traits.keys():
        taxalist.append(int(s))

    inherited = []

    results = {}
    n_mutations_d = []
    n_mutations_c = []

    all_focal_trees = []
    counts_by_tree = []

    
    events = []
    for e in admix:
        events.append([e[0], str(conversions[e[1]]), str(conversions[e[2]]), e[3]])
    admix = events

    if len(admix) != 0:
        total_reps_for_intro = 0
        for e in admix:
            total_reps_for_intro += int(reps * float(e[3]))
    remaining_reps = reps - total_reps_for_intro


    per_thread = remaining_reps//threads
    v = remaining_reps/threads
    per_thread = [per_thread]*threads
    if not v.is_integer():
        threads += 1
        per_thread.append(remaining_reps%(threads-1))


    if len(admix) != 0:
        #Extra thread for introgression
        threads += 1




    prefix = args.outputdir

    processes_ms = []
    processes_sq = []
    intro_indices = []
    if len(admix) == 0:
        for y in range(0, threads):
            ms_call = hemiplasytool.splits_to_ms(splits, taxa, per_thread[y], args.mspath, y, prefix)
            m = hemiplasytool.call_programs(ms_call, "", "trees.tmp", taxalist)
            processes_ms.append(m)
    elif len(admix) != 0:
        for y in range(0, threads):
            if y != threads-1:
                ms_call = hemiplasytool.splits_to_ms(splits, taxa, per_thread[y], args.mspath, y, prefix)
                m = hemiplasytool.call_programs(ms_call, "", "trees.tmp", taxalist)
                processes_ms.append(m)
            elif y == threads-1:
                ms_calls = []
                for m, event in enumerate(admix):
                    o = str(y) + "_" + str(m)
                    intro_indices.append(m)
                    ms_call = hemiplasytool.splits_to_ms(splits, taxa, 
                        int(reps * float(event[3])), args.mspath, o, prefix, event)
                    m = hemiplasytool.call_programs(ms_call, "", "trees.tmp", taxalist)
                    processes_ms.append(m)

    done = False
    while done == False:
        ms_processes = []
        for p in processes_ms:
            poll = p.poll()
            if poll == None:
                ms_processes.append(False)
            else:
                ms_processes.append(True)
    
        j = all(process == True for process in ms_processes)
        done = j



    string_cat_ms = "cat "
    for y in range(0, threads):
        if y != threads-1:
            string_cat_ms += prefix + ".trees" + str(y) + ".tmp "
        elif y == threads-1:
            for intro in intro_indices:
                string_cat_ms += prefix + ".trees" + str(y) + "_" + str(intro) + ".tmp "
    string_cat_ms += "> " + prefix + ".trees.tmp"
    os.system(string_cat_ms)

    for y in range(0, threads):
        if y != threads-1:
            seqgencall = hemiplasytool.seq_gen_call(prefix + ".trees" + str(y) + ".tmp", args.seqgenpath, args.mutationrate, str(y), prefix)
            s = hemiplasytool.call_programs_sg(ms_call, seqgencall, "trees.tmp", taxalist)
            processes_sq.append(s)
        else:
            for intro in intro_indices:
                seqgencall = hemiplasytool.seq_gen_call(prefix + ".trees" + str(y) + "_" + str(intro) + ".tmp", args.seqgenpath, args.mutationrate, str(y), prefix)
                s = hemiplasytool.call_programs_sg(ms_call, seqgencall, "trees.tmp", taxalist)
                processes_sq.append(s)


    intro_start = sum(per_thread)

    done = False
    while done == False:
        sg_processes = []
        for p in processes_sq:
            poll = p.poll()
            if poll == None:
                sg_processes.append(False)
            else:
                sg_processes.append(True)
        j = all(process == True for process in sg_processes)
        done = j

    string_cat = "cat "
    for y in range(0, threads):
        string_cat += prefix + ".seqs" + str(y) + ".tmp "
    string_cat += "> " + prefix + ".seqs.tmp"
    os.system(string_cat)

    # Gets indices of trees with site patterns that match speecies pattern
    log.debug("Finding trees that match species trait pattern...")
    match_species_pattern, counts_by_tree = seqtools.readSeqs(
        prefix + ".seqs.tmp", len(taxalist), traits, len(splits), i, prefix, intro_start
    )

    log.debug("Getting focal trees...")
    # Gets the trees at these indices 
    focal_trees, _ = seqtools.getTrees(prefix + ".trees.tmp", match_species_pattern)
    all_focal_trees = focal_trees
    assert len(match_species_pattern) == len(focal_trees)
    log.debug("Calculating discordance...")
    results[i], disc, conc = seqtools.propDiscordant(focal_trees, treeSp)

    focaltrees_d = seqtools.parse_seqgen(prefix + ".focaltrees.tmp", len(taxalist), disc)
    focaltrees_c = seqtools.parse_seqgen(prefix + ".focaltrees.tmp", len(taxalist), conc)
    for index, tree in enumerate(focaltrees_d):
        n_mutations_d.append(seqtools.count_mutations(tree, len(taxalist)))
    for index, tree in enumerate(focaltrees_c):
        n_mutations_c.append(seqtools.count_mutations(tree, len(taxalist)))
    nderived = 0
    for trait in traits.values():
        if trait == 1:
            nderived += 1
    interesting = seqtools.get_interesting(
        focaltrees_d, nderived, len(traits.keys())
    )
    for item in interesting:
        test_summarize = seqtools.summarize_interesting(item, len(traits.keys()))
        inherited = inherited + test_summarize

    # Clean up temporary files
    os.system("rm *.tmp")
    ###################################################################

    # Begin summary of all batches
    mutation_counts_d = [[x, n_mutations_d.count(x)] for x in set(n_mutations_d)]
    mutation_counts_c = [[x, n_mutations_c.count(x)] for x in set(n_mutations_c)]
    summary = hemiplasytool.summarize(results)
    #counts_by_tree = seqtools.sum_counts_by_tree(counts_by_tree)
    if len(inherited) > 0:
        mutation_pat = hemiplasytool.summarize_inherited(inherited)
    else:
        mutation_pat = None
        log.debug(
            "Not enough 'interesting' cases to provide mutation inheritance patterns"
        )
    min_mutations_required = hemiplasytool.fitchs_alg(str(treeSp), traits)

    #log.debug("Plotting...")
    #try:
    #    hemiplasytool.plot_mutations(mutation_counts_c, mutation_counts_d, args.outputdir)
    #except:
    #    log.debug("Can't plot!")


    log.debug("Writing output file...")
    hemiplasytool.write_output(
        summary,
        mutation_counts_c,
        mutation_counts_d,
        mutation_pat,
        counts_by_tree,
        str(treeSp),
        admix,
        traits,
        min_mutations_required,
        args.outputdir,
        (reps),
        conversions,
        original_tree[0],
        intercept,
        coef,
        newick_internals,
        coal_internals
    )
    hemiplasytool.write_unique_trees(all_focal_trees, args.outputdir, traits)
    end = time.time()
    print("\nTime elapsed: " + str(end - start) + " seconds")
    ################################################################

if __name__ == "__main__":
    main(*sys.argv)
