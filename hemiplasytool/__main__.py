# /usr/bin/python3
import argparse
import time
import sys
import logging as log
from hemiplasytool import hemiplasytool
from hemiplasytool import seqtools

"""
HemiplasyTool
Authors: Matt Gibson, Mark Hibbins
Indiana University
"""


def main(*args):
    start = time.time()
    parser = argparse.ArgumentParser(description="Tool for characterising \
                        hemiplasy given traits mapped onto a species tree")
    parser.add_argument("-v", "--verbose",
                        help="Enable debugging messages to be displayed",
                        action='store_true')
    parser.add_argument("input", metavar='splits',
                        help="Input file describing split times, trait pattern, and topology")
    parser.add_argument("-n", "--replicates", metavar="",
                        help="Number of replicates per batch", default=1000000)
    parser.add_argument("-x", "--batches", metavar="",
                        help="Number of batches", default=3)
    parser.add_argument("-p", "--mspath", metavar="",
                        help="Path to ms", default="./msdir")
    parser.add_argument("-g", "--seqgenpath", metavar="",
                        help="Path to seq-gen", default="./seq-gen")
    parser.add_argument("-s", "--mutationrate", metavar="",
                        help="Seq-gen mutation rate (default 0.05)",
                        default=0.05)
    parser.add_argument("-o", "--outputdir", metavar="",
                        help="Output directory")

    args = parser.parse_args()

    # Setup logging
    log.basicConfig(level=log.DEBUG)
    logger = log.getLogger()
    if (args.verbose):
        logger.disabled = False
    else:
        logger.disabled = True

    # Silence matplotlib debug logging
    mpl_logger = log.getLogger('matplotlib')
    mpl_logger.setLevel(log.WARNING)

    # Read input file
    log.debug("Reading input file...")
    splits, taxa, traits, speciesTree,\
        admix = hemiplasytool.readInput(args.input)

    batches = int(args.batches)
    reps = int(args.replicates)

    # Make program calls
    if len(admix) != 0:
        log.debug("Introgression events specified")
        ms_call = []
        summ = 0.0

        # Intogression trees
        for i, event in enumerate(admix):
            summ += float(event[3])
            ms_call.append(hemiplasytool.splits_to_ms(splits, taxa, int(reps*float(event[3])), args.mspath, event, i))
        
        # Species tree
        ms_call.append(hemiplasytool.splits_to_ms(splits, taxa, int(reps*(1-summ)), args.mspath, event, i+1))
    else:
        ms_call = hemiplasytool.splits_to_ms(splits, taxa, args.replicates, args.mspath)

    seqgencall = hemiplasytool.seq_gen_call('trees.tmp', args.seqgenpath, args.mutationrate)

    taxalist = []
    for s in traits.keys():
        taxalist.append(int(s))

    inherited = []

    results = {}
    n_mutations_d = []
    n_mutations_c = []
    for i in range(0, batches):


        hemiplasytool.call_programs(ms_call, seqgencall, 'trees.tmp', taxalist)

        # Gets indices of trees with site patterns that match speecies pattern
        log.debug("Extracting trees that match species trait pattern...")
        match_species_pattern = seqtools.readSeqs("seqs.tmp", len(taxalist), traits, len(splits), i)

        # Gets the trees at these indices
        focal_trees, _ = seqtools.getTrees('trees.tmp', match_species_pattern)

        assert len(match_species_pattern) == len(focal_trees)

        """
        Out of those trees which follow the species site pattern,
        get the number of trees which are discordant.
        """

        log.debug("Calculating discordance...")
        results[i], disc, conc = seqtools.propDiscordant(focal_trees, speciesTree)

        focaltrees_d = seqtools.parse_seqgen("focaltrees.tmp", 
            len(taxalist), disc)
        focaltrees_c = seqtools.parse_seqgen("focaltrees.tmp", 
            len(taxalist), conc)

        for index, tree in enumerate(focaltrees_d):
            n_mutations_d.append(seqtools.count_mutations(tree, len(taxalist)))
        for index, tree in enumerate(focaltrees_c):
            n_mutations_c.append(seqtools.count_mutations(tree, len(taxalist)))

        nderived = 0
        for trait in traits.values():
            if trait == '1':
                nderived += 1

        interesting = seqtools.get_interesting(focaltrees_d, nderived, len(traits.keys()))
        for item in interesting:
            test_summarize = seqtools.summarize_interesting(item, len(traits.keys()))
            inherited = inherited + test_summarize

        # Clean up temporary files from this batch
        hemiplasytool.cleanup()

    mutation_counts_d = [[x,n_mutations_d.count(x)] for x in set(n_mutations_d)]
    mutation_counts_c = [[x,n_mutations_c.count(x)] for x in set(n_mutations_c)]

    summary = hemiplasytool.summarize(results)

    print("\n####################RESULTS####################\n###############################################")

    print("\nOf the replicates that follow species site pattern: ")
    print(str(summary[0]) + " were discordant\n" + str(summary[1]-summary[0]) + " were concordant\n")

    print("\nOn concordant trees:")
    print("# Mutations\t# Trees")
    for item in mutation_counts_c:
        print(str(item[0]) + '\t\t' + str(item[1]))
    print("\nOn discordant trees:")
    print("# Mutations\t# Trees")
    for item in mutation_counts_d:
        print(str(item[0]) + '\t\t' + str(item[1]))

    if len(inherited) > 0:
        mutation_pat = hemiplasytool.summarize_inherited(inherited)
    else:
        mutation_pat = None
        log.debug("Not enough 'interesting' cases to provide mutation inheritance patterns")

    print()
    log.debug("Plotting...")
    hemiplasytool.plot_mutations(mutation_counts_c, mutation_counts_d)

    hemiplasytool.write_output(summary, mutation_counts_c, mutation_counts_d, mutation_pat, args.outputdir)

    end = time.time()
    print("\nTime elapsed: " + str(end - start) + " seconds")


if __name__ == "__main__":
    main(*sys.argv)
