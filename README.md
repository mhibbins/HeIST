# HemiplasyTool

## Authors:
Matt Gibson (gibsomat@indiana.edu)  
Mark Hibbins (mhibbins@indiana.edu)

## Dependencies:
* ms  
* seq-gen  
* biopython
* numpy


## Installation

### Install from source
```
git clone https://github.com/mhibbins/hemiplasytool
cd hemiplasytool
python setup.py install
```

### Install from PyPI
```
pip install hemiplasytool
```

## Usage
```
usage: hemiplasytool [-h] [-v] [-n] [-x] [-p] [-g] [-s] [-o] splits

Calculate the probability that convergent trait patterns are due to hemiplasy

positional arguments:
  splits                Input file describing split times, trait pattern, and
                        topology

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Enable debugging messages to be displayed
  -n , --replicates     Number of replicates per batch
  -x , --batches        Number of batches
  -p , --mspath         Path to ms
  -g , --seqgenpath     Path to seq-gen
  -s , --mutationrate   Seq-gen mutation rate (default 0.05)
  -o , --outputdir      Output directory
```

## Input file

The input file has three sections:  split times, traits, and species tree. They must be specified in this order and delimited by a '#'. See below for descriptions of each section

```
#splits
6   2   1
3   3   2
1.5 5   3
1.25    6   5
1   4   3

#traits
1   0 
2   1
3   0
4   1
5   0
6   1

#tree
(1,(2,((6,5),(4,3))));

#admixture (time, source, dest, probability)
0.25    3   2   0.1
0.5 5   6   0.1

```

### Split times

The split times describe the order of subpopulation splits to `ms`. Each line specifies the timing (in 4N generations), source population, and destination population (backwards in time). Splits should be ordered oldest to newest. Entries should be delimited by spaces or tabs


### Traits

The traits section describes the observed species trait pattern. Each line specifies the taxa ID (must correspond to those coded in the split times file), the binary trait value, and the timing of sampling (in 4N generations relative to the longest branch). These can be specified in any order


### Species tree

The species tree in Newick format. Again, taxa IDs must correspond to those in the split times and traits sections.


## Example:
```
hemiplasytool -v -n 1000000 -p ~/bin/ms -g ~/bin/seq-gen -x 1 ./input_test.txt
```

### Output:
```
Of the replicates that follow species site pattern:
118 were discordant
32 were concordant


On concordant trees:
# Mutations	# Trees
3		28
4		3
5		1

On discordant trees:
# Mutations	# Trees
1		5
2		21
3		70
4		20
5		2

Derived mutation inheritance patterns for trees with fewer mutations than derived taxa:

	Term	Inherited from anc node
Taxa 2	15	6
Taxa 4	0	21
Taxa 6	1	20

DEBUG:root:Plotting...

Time elapsed: 47.09378099441528 seconds
```

![Mutation distribution](mutation_dist.png)