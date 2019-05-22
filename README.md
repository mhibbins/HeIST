# HemiplasyTool


![](https://img.shields.io/github/release-pre/mhibbins/hemiplasytool.svg)
![](https://img.shields.io/github/release-date-pre/mhibbins/hemiplasytool.svg)
## Authors:
Matt Gibson (gibsomat@indiana.edu)  
Mark Hibbins (mhibbins@indiana.edu)

## Dependencies:
* [ms](http://home.uchicago.edu/~rhudson1/source.html)  
* [seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/)
* biopython
* numpy
* matplotlib


## Installation

```
git clone https://github.com/mhibbins/hemiplasytool
cd hemiplasytool
python setup.py install
```


## Usage
```
usage: hemiplasytool [-h] [-v] [-n] [-x] [-p] [-g] [-s] [-o] splits

Tool for characterising hemiplasy given traits mapped onto a species tree

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

#introgression (time, source, dest, probability; optional)
0.25    3   2   0.1
0.5 5   6   0.1

```

### Split times

The split times describe the order of subpopulation splits to `ms`. Each line specifies the timing (in 4N generations), source population, and destination population (backwards in time). Splits should be ordered oldest to newest. Entries should be delimited by spaces or tabs


### Traits

The traits section describes the observed species trait pattern. Each line specifies the taxa ID (must correspond to those coded in the split times file), the binary trait value, and the timing of sampling (in 4N generations relative to the longest branch). These can be specified in any order


### Species tree

The species tree in Newick format. Again, taxa IDs must correspond to those in the split times and traits sections.

### Introgression

Introgression events. Each line should specify the timing (in 4N generations), source taxon, destination taxon, and probability of introgression. Events can be specified in any order.

## Example:
```
python -m hemiplasytool -v -n 1000000 -p ~/bin/ms -g ~/bin/seq-gen -x 1 ./test/input_test_small.txt -o outtest.txt
```

### Output:
```
OOn concordant trees:
# Mutations	# Trees
3		25
4		2
5		2

On discordant trees:
# Mutations	# Trees
1		5
2		23
3		79
4		9
5		3

Derived mutation inheritance patterns for trees with fewer         mutations than derived taxa:

	Term	Inherited from anc node
Taxa 2	9	14
Taxa 4	1	22
Taxa 6	2	21

Of the replicates that follow species site pattern:
119 were discordant
29 were concordant

9 replicates matching the species pattern were from introgression tree 1
18 replicates matching the species pattern were from introgression tree 2
121 replicates matching the species pattern were from the species tree

DEBUG:root:Plotting...
DEBUG:root:Writing output file...

Time elapsed: 50.345940828323364 seconds
```

See `outtest.txt` for full output.

![Mutation distribution](mutation_dist.png)