# HemiplasyTool

## Dependencies:
* ms  
* seq-gen  
* ete3   

### Install ete3 with Anaconda
```
conda install -c etetoolkit ete3 ete_toolchain
```

## Usage
```
usage: hemiplasytool.py [-h] [-v] [-s SPLITTIMES] [-t TRAITS] [-b SPECIESTREE]
                        [-n REPLICATES] [-x BATCHES] [-p MSPATH]
                        [-g SEQGENPATH] [-o OUTPUTDIR]

Calculate the probability that convergent trait patterns are due to hemiplasy

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Enable debugging messages to be displayed
  -s SPLITTIMES, --splittimes SPLITTIMES
                        Split times file, ordered from oldest to newest. In
                        units of 4N generations.
  -t TRAITS, --traits TRAITS
                        Traits file
  -b SPECIESTREE, --speciestree SPECIESTREE
                        Species topology in Newick format on one line.
  -n REPLICATES, --replicates REPLICATES
                        Number of replicates per batch
  -x BATCHES, --batches BATCHES
                        Number of batches
  -p MSPATH, --mspath MSPATH
                        Path to ms
  -g SEQGENPATH, --seqgenpath SEQGENPATH
                        Path to seq-gen
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        Output directory
```

## Input files

### Split times file

The split times file describes the order of subpopulation splits to `ms`. Each line in the file specifies the timing (in 4N generations), source population, and destination population (backwards in time). Splits should be ordered oldest to newest. Entries should be delimited by spaces or tabs

```
6   2   1
3   3   2
1.5 5   3
1.25    6   5
1   4   3
```

### Traits file

The traits file describes the observed species trait pattern. Each line in the file specifies the taxa ID (must correspond to those coded in the split times file), the binary trait value, and the timing of sampling (in 4N generations relative to the longest branch).

```
1   0   0.3
2   1   0.2
3   0   0.05
4   1   0
5   0   0.1
6   1   0.15
```

### Species tree file

The species tree file specifies the tree topology in Newick format. Again, taxa IDs must correspond to those in the split times and traits file.

```
(1,(2,((6,5),(4,3))));
```


## Example:
```
python hemiplasytool.py -s ../splits_test.txt -n 1000000 -x 2 -p ~/bin/ms -t ../traits_test.txt -g  ~/bin/seq-gen -b ../test_species.txt -v
```