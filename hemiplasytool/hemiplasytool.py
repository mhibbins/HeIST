"""
Hemiplasy Tool
Authors: Matt Gibson, Mark Hibbins
Indiana University
"""
import argparse
import sys
import os
import logging as log
import seqtools
import time
import matplotlib.pyplot as plt
import numpy as np


def splits_to_ms(splitTimes, taxa, reps, path_to_ms):
    """
    Converts inputs into a call to ms
    TODO: Add introgression
    """
    nsamples = len(splitTimes)+1
    call = path_to_ms + ' ' + str(nsamples) + ' ' + str(reps) + ' -T -I ' + str(nsamples) + ' '
    for i in range(0, nsamples):
        call += '1 '
    for x, split in enumerate(splitTimes):
        call += '-ej ' + str(split) + ' ' + str(taxa[x][0]) + ' ' + str(taxa[x][1]) + ' '

    call += "| tail -n +4 | grep -v // > trees.tmp"
    return(call)

def seq_gen_call(treefile, path, s=0.05):
    """
    Make seq-gen call.
    """
    return(path + ' -m HKY -l 1 -s ' + str(s) + ' -wa <"' + treefile + '" > seqs.tmp')

def call_programs(ms_call, seqgencall, treefile, ntaxa):
    """
    Calls ms and seq-gen
    """
    log.debug("Calling ms...")
    os.system(ms_call)

    #sedTrees(treefile, ntaxa)

    log.debug("Calling seq-gen...")
    os.system(seqgencall)
    #os.system("rm trees.tmp; rm trees.tmp.bak")

def sedTrees(treefile, taxalist):
    """
    Flips taxa IDs with sed. They become reversed when taking "ancient" samples.
    Currently not used.**
    """
    newTaxa = {}
    max1 = max(taxalist)
    for t in taxalist:
        newTaxa[t] = max1 - t + 1
    seen = []
    for key,val in newTaxa.items():
        if (key not in seen):
            seen.append(key)
            seen.append(val)
            if sys.platform == "darwin":
                call = "sed -i '.bak' 's/[[:<:]]" + str(key) + ":[[:>:]]/~~" + "/g; s/[[:<:]]" + str(val) + ":[[:>:]]/[[:<:]]" + str(key) + ":[[:>:]]/g; s/~~/[[:<:]]" + str(val) + ":[[:>:]]/g' " + treefile 
                #print(call)
                log.debug("Fixing taxa names...")
                os.system(call)
            elif sys.platform == "linux" or sys.platform == "linux2":
                call = "sed -i 's/" + str(key) + "/~~" + "/g; s/" + str(val) + "/" + str(key) + "/g; s/~~/" + str(val) + "/g' " + treefile 
                log.debug("Fixing taxa names...")
                os.system(call)
    if sys.platform == "darwin":
        os.system("rm " + treefile + ".bak")

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

def write_output(summary, mutation_counts_c, mutation_counts_d, filename):
    out1 = open(filename, 'w')

    out1.write("Of the replicates that follow species site pattern:\n")
    out1.write(str(summary[0]) + " were discordant\n" + str(summary[1]-summary[0]) + " were concordant\n")

    out1.write("\nOn concordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in  mutation_counts_c:
        out1.write(str(item[0]) + '\t\t' + str(item[1]) + '\n')
    out1.write("\nOn discordant trees:\n")
    out1.write("# Mutations\t# Trees\n")
    for item in  mutation_counts_d:
        out1.write(str(item[0]) + '\t\t' + str(item[1]) + '\n')
    out1.close()

def plot_mutations(results_c, results_d):
    """
    Plot mutation distribution with matplotlib
    """
    objs_c = [i[0] for i in results_c]
    objs_d = [i[0] for i in results_d]
    conc_dic = {}
    disc_dic = {}
    objs = objs_c  + objs_d
    objs = set(objs)
    x = np.array(list(range(1,max(objs)+1)))
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
    
#def fishers_exact(counts):
#    """Scipy fishers exact test"""
#    return(fisher_exact([[counts[0], (counts[1]-counts[0])],[counts[2], (counts[3]-counts[2])]]))


def readInput(file):
    f = open(file, 'r')
    splits = []
    taxa = []
    traits = {}

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
            tree = line.replace("\n","")
    return(splits, taxa, traits, tree)

    
def main(*args):
    start = time.time()
    parser = argparse.ArgumentParser(description="Calculate the probability that convergent trait patterns are due to hemiplasy")
    parser.add_argument("-v", "--verbose", help="Enable debugging messages to be displayed", action='store_true')
    parser.add_argument("input", metavar='splits', help="Input file describing split times, trait pattern, and topology")
    parser.add_argument("-n","--replicates", metavar="", help="Number of replicates per batch", default=1000000)
    parser.add_argument("-x","--batches", metavar="", help="Number of batches", default=3)
    parser.add_argument("-p","--mspath", metavar="", help="Path to ms", default="./msdir")
    parser.add_argument("-g","--seqgenpath", metavar="", help="Path to seq-gen", default="./seq-gen")
    parser.add_argument("-s","--mutationrate", metavar="", help="Seq-gen mutation rate (default 0.05)", default=0.05)
    parser.add_argument("-o","--outputdir", metavar="", help="Output directory")
    args = parser.parse_args()

    #Setup logging
    log.basicConfig(level=log.DEBUG)
    logger = log.getLogger()
    if (args.verbose):
        logger.disabled = False
    else:
        logger.disabled =  True
    #Silence matplotlib debug logging
    mpl_logger = log.getLogger('matplotlib')
    mpl_logger.setLevel(log.WARNING) 


    #read input file
    log.debug("Reading input file...")
    splits, taxa, traits, speciesTree = readInput(args.input)

    batches = int(args.batches)

    #Make program calls
    ms_call = splits_to_ms(splits, taxa, args.replicates, args.mspath)
    seqgencall = seq_gen_call('trees.tmp', args.seqgenpath, args.mutationrate)

    log.debug(ms_call)
    log.debug(seqgencall)


    taxalist = []
    for s in traits.keys():
        taxalist.append(int(s))

    results = {}
    results_alltrees ={}
    n_mutations_d = []
    n_mutations_c = []
    for i in range(0, batches):
        #Call ms and seq-gen
        call_programs(ms_call, seqgencall, 'trees.tmp', taxalist)

        #Gets indices of trees with site patterns that match speecies pattern
        log.debug("Extracting trees that match species trait pattern...")
        match_species_pattern = seqtools.readSeqs("seqs.tmp",len(taxalist), traits, len(splits), i)

        #Gets the trees at these indices
        focal_trees, _ = seqtools.getTrees('trees.tmp', match_species_pattern)

        assert len(match_species_pattern) == len(focal_trees)
        out3 = open('focal_newick.tmp', 'w')
        for z in focal_trees:
            out3.write(z + '\n')
        out3.close()

        #Out of those trees which follow the species site pattern, get the number
        #of trees which are discordant.
        log.debug("Calculating discordance...")
        results[i], disc, conc = seqtools.propDiscordant(focal_trees, speciesTree)
        #TODO: Add catch here. If # that follow is very low, restart loop with higher value for n

        #log.debug("Calculating discordance...")
        #results_alltrees[i], _, _ = seqtools.propDiscordant_async(all_trees, speciesTree)
        
        focaltrees_d = seqtools.parse_seqgen("focaltrees.tmp", len(taxalist), disc)
        focaltrees_c = seqtools.parse_seqgen("focaltrees.tmp", len(taxalist), conc)

        for index, tree in enumerate(focaltrees_d):
	        n_mutations_d.append(seqtools.count_mutations(tree, len(taxalist)))
        for index, tree in enumerate(focaltrees_c):
            n_mutations_c.append(seqtools.count_mutations(tree, len(taxalist)))
        
        #Clean up temporary files from this batch
        cleanup()
        
    mutation_counts_d = [[x,n_mutations_d.count(x)] for x in set(n_mutations_d)]
    mutation_counts_c = [[x,n_mutations_c.count(x)] for x in set(n_mutations_c)]

    summary = summarize(results)

    print("\n####################RESULTS####################\n###############################################")

    print("\nOf the replicates that follow species site pattern: ")
    print(str(summary[0]) + " were discordant\n" + str(summary[1]-summary[0]) + " were concordant\n")


  
    print("\nOn concordant trees:")
    print("# Mutations\t# Trees")
    for item in  mutation_counts_c:
        print(str(item[0]) + '\t\t' + str(item[1]))
    print("\nOn discordant trees:")
    print("# Mutations\t# Trees")
    for item in  mutation_counts_d:
        print(str(item[0]) + '\t\t' + str(item[1]))
    
    log.debug("Plotting...")
    plot_mutations(mutation_counts_c, mutation_counts_d)


    write_output(summary, mutation_counts_c, mutation_counts_d, args.outputdir)

    end = time.time()
    print("\nTime elapsed: " + str(end - start) + " seconds")

if __name__ == "__main__":
    main(*sys.argv)

