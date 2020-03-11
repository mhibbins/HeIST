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
from hemiplasytool import hemiplasytool
from hemiplasytool import seqtools
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

def main(*args):
    start = time.time()
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


    # Convert coalescent tree to ms splits
    treeSp, conversions = hemiplasytool.names2ints(treeSp)
    original_tree[0], _ = hemiplasytool.names2ints(original_tree[0])
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

    breaks = []

    # Convert introgression taxa to ints
    events = []
    for e in admix:
        events.append([e[0], str(conversions[e[1]]), str(conversions[e[2]]), e[3]])
    admix = events
    if len(admix) != 0:
        breaks = [0] * len(admix)
        log.debug("Introgression events specified")
        ms_call = []
        summ = 0.0
        # Intogression trees
        for i, event in enumerate(admix):
            summ += float(event[3])
            if i == 0:
                breaks[i] = int(reps * float(event[3]))
            else:
                breaks[i] = int(reps * float(event[3])) + breaks[i - 1]
            ms_call.append(
                hemiplasytool.splits_to_ms(
                    splits, taxa, int(reps * float(event[3])), args.mspath, event, i
                )
            )

        # Species tree
        ms_call.append(
            hemiplasytool.splits_to_ms(
                splits, taxa, int(reps * (1 - summ)), args.mspath, event, i + 1
            )
        )
        breaks.append(int(reps * (1 - summ)) + breaks[len(breaks) - 1])
    else:
        pass
        #ms_call = hemiplasytool.splits_to_ms(splits, taxa, args.replicates, args.mspath)

    #seqgencall = hemiplasytool.seq_gen_call(
    #    "trees.tmp", args.seqgenpath, args.mutationrate
    #)
    #################################################################


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

    per_thread = reps//threads
    v = reps/threads
    per_thread = [per_thread]*threads
    if not v.is_integer():
        threads += 1
        per_thread.append(reps%(threads-1))

    processes_ms = []
    processes_sq = []
    for y in range(0, threads):
        ms_call = hemiplasytool.splits_to_ms(splits, taxa, per_thread[y], args.mspath, y)
        m = hemiplasytool.call_programs(ms_call, "", "trees.tmp", taxalist)
        processes_ms.append(m)

    print("RUNNING MS!")
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

    print("CONCAT MS!")
    string_cat_ms = "cat "
    for y in range(0, threads):
        string_cat_ms += "trees" + str(y) + ".tmp "
    string_cat_ms += "> trees.tmp"
    os.system(string_cat_ms)

    print("RUNNING SEQ GEN!")
    for y in range(0, threads):
        seqgencall = hemiplasytool.seq_gen_call("trees" + str(y) + ".tmp", args.seqgenpath, args.mutationrate, str(y))
        #print(seqgencall)
        s = hemiplasytool.call_programs_sg(ms_call, seqgencall, "trees.tmp", taxalist)
        processes_sq.append(s)

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

    print("CONCAT SEQ GEN!")
    string_cat = "cat "
    for y in range(0, threads):
        string_cat += "seqs" + str(y) + ".tmp "
    string_cat += "> seqs.tmp"
    os.system(string_cat)
    print(str(time.time() - start))

    print("READSEQS!")
    # Gets indices of trees with site patterns that match speecies pattern
    log.debug("Extracting trees that match species trait pattern...")
    match_species_pattern, counts = seqtools.readSeqs(
        "seqs.tmp", len(taxalist), traits, len(splits), i, breaks
    )
    counts_by_tree.append(counts)
    print(str(time.time() - start))
    print("GET TREES!")
    # Gets the trees at these indices
    focal_trees, _ = seqtools.getTrees("trees.tmp", match_species_pattern)
    all_focal_trees = all_focal_trees + focal_trees
    print(str(time.time() - start))

    assert len(match_species_pattern) == len(focal_trees)

    log.debug("Calculating discordance...")
    print("PROP DISCORDANT!")
    results[i], disc, conc = seqtools.propDiscordant(focal_trees, treeSp)
    print(str(time.time() - start))

    print("PARSE SEQGEN!")
    focaltrees_d = seqtools.parse_seqgen("focaltrees.tmp", len(taxalist), disc)
    focaltrees_c = seqtools.parse_seqgen("focaltrees.tmp", len(taxalist), conc)
    print(str(time.time() - start))

    print("COUNT MUTATIONS!")
    for index, tree in enumerate(focaltrees_d):
        n_mutations_d.append(seqtools.count_mutations(tree, len(taxalist)))
    for index, tree in enumerate(focaltrees_c):
        n_mutations_c.append(seqtools.count_mutations(tree, len(taxalist)))
    print(str(time.time() - start))

    nderived = 0
    for trait in traits.values():
        if trait == 1:
            nderived += 1
    print("GET INTERESTING!")
    interesting = seqtools.get_interesting(
        focaltrees_d, nderived, len(traits.keys())
    )
    print(str(time.time() - start))

    print("SUMMARIZE INTERESTING!!")
    for item in interesting:
        test_summarize = seqtools.summarize_interesting(item, len(traits.keys()))
        inherited = inherited + test_summarize
    print(str(time.time() - start))

    # Clean up temporary files
    os.system("rm *.tmp")
    ###################################################################


    # Begin summary of all batches
    mutation_counts_d = [[x, n_mutations_d.count(x)] for x in set(n_mutations_d)]
    mutation_counts_c = [[x, n_mutations_c.count(x)] for x in set(n_mutations_c)]

    summary = hemiplasytool.summarize(results)

    counts_by_tree = seqtools.sum_counts_by_tree(counts_by_tree)

    if len(inherited) > 0:
        mutation_pat = hemiplasytool.summarize_inherited(inherited)
    else:
        mutation_pat = None
        log.debug(
            "Not enough 'interesting' cases to provide mutation inheritance patterns"
        )

    min_mutations_required = hemiplasytool.fitchs_alg(str(treeSp), traits)

    log.debug("Plotting...")
    try:
        hemiplasytool.plot_mutations(mutation_counts_c, mutation_counts_d, args.outputdir)
    except:
        log.debug("Can't plot!")


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
