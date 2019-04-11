"""
Hemiplasy Tool
Authors: Matt Gibson, Mark Hibbins
Indiana University
"""
import argparse
import sys
import os
import logging as log

def read_splits(file):
    """
    Splits file should have three values per line: Timing of split (4N gens) \t Source pop \t Dest pop.
    Splits should be ordered oldest to most recent.
    """
    input1 = open(file, 'r')
    splits = []
    taxa = []
    for i, line in enumerate(input1):
        l = line.replace("\n", "").split()
        splits.append(float(l[0]))
        taxa.append((int(l[1]), int(l[2])))
    input1.close()
    return(splits, taxa)

def read_traits(file):
    """
    Traits file should have three values per line: Taxa # \t Binary trait value (0 or 1) \t Timing of sampling, if tree is non-ultrametric, relative to longest terminal branch. 
    The longest terminal branch should have a value of 0. If species tree is ultrametric, please specify all 0s in the third column.
    """
    input1 = open(file, 'r')
    traits = {}
    sample_times = {}
    for i, line in enumerate(input1):
        l = line.replace('\n','').split()
        traits[l[0]] = l[1]
        sample_times[l[0]] = l[2]
    return(traits, sample_times)

def splits_to_ms(splitTimes, taxa, reps, sampleTimes, path_to_ms):
    """
    Converts inputs into a call to ms
    """
    for key, val in sampleTimes.items():
        if (val == '0') or (val == '0.0'):
            index_of_most_recent = int(key)

    nsamples = len(splitTimes)+1
    call = path_to_ms + ' ' + '1' + ' ' + str(reps) + ' -T -I ' + str(nsamples) + ' '
    for i in range(0, nsamples):
        if (i+1) == index_of_most_recent:
            call += '1 '
        else:
            call += '0 '
    for x, split in enumerate(splitTimes):
        call += '-ej ' + str(split) + ' ' + str(taxa[x][0]) + ' ' + str(taxa[x][1]) + ' '
    
    for time in sampleTimes.values():
        if time != '0':
            for taxa, time in sampleTimes.items():
                if (float(time) > 0):
                    call += '-eA ' + time + ' ' + taxa + ' 1 '
        break
    call += "| tail -n +4 | grep -v // > trees.tmp"
    return(call)

def seq_gen_call(treefile, path):
    """
    Make seq-gen call. TODO: Add option to change seq-gen parameters.
    """
    return(path + ' -m HKY -l 1 -s 0.01 -wa <"' + treefile + '" > seqs.tmp')

def call_programs(ms_call, seqgencall, treefile, ntaxa):
    """
    Calls ms and seq-gen
    """
    log.debug("Calling ms...")
    os.system(ms_call)

    sedTrees(treefile, ntaxa)

    log.debug("Calling seq-gen...")
    os.system(seqgencall)
    #os.system("rm trees.tmp; rm trees.tmp.bak")

def sedTrees(treefile, taxalist):
    """
    Flips taxa IDs with sed. They become reversed when taking "ancient" samples.
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
	    #TODO: deal with .bak call for mac vs. linux systems
            call = "sed -i '.bak' 's/" + str(key) + "/~~" + "/g; s/" + str(val) + "/" + str(key) + "/g; s/~~/" + str(val) + "/g' " + treefile 
            log.debug("Fixing taxa names...")
            os.system(call)

def main(*args):
    parser = argparse.ArgumentParser(description="Calculate the probability that convergent trait patterns are due to hemiplasy")
    parser.add_argument("-v", "--verbose", help="Enable debugging messages to be displayed", action='store_true')
    parser.add_argument("-s","--splittimes", help="Split times file, ordered from oldest to newest. In units of 4N generations.")
    parser.add_argument("-t","--traits", help="Traits file")  
    parser.add_argument("-n","--replicates", help="Number of replicates")
    parser.add_argument("-p","--mspath", help="Path to ms")
    parser.add_argument("-g","--seqgenpath", help="Path to seq-gen")
    parser.add_argument("-o","--outputdir", help="Output directory")
    args = parser.parse_args()

    #Setup logging
    log.basicConfig(level=log.DEBUG)
    logger = log.getLogger()
    if (args.verbose):
        logger.disabled = False
    else:
        logger.disabled =  True
    
    #read input files
    splits, taxa = read_splits(args.splittimes)
    traits, sample_times = read_traits(args.traits)

    #Make program calls
    ms_call = splits_to_ms(splits, taxa, args.replicates, sample_times, args.mspath)
    seqgencall = seq_gen_call('trees.tmp', args.seqgenpath)

    taxalist = []
    for s in sample_times.keys():
        taxalist.append(int(s))
    
    #Call ms and seq-gen
    call_programs(ms_call, seqgencall, 'trees.tmp', taxalist)



if __name__ == "__main__":
    main(*sys.argv)

