# HeIST
**He**miplasy **I**nference **S**imulation **T**ool


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
* ete3


## Installation

### From GitHub
```
git clone https://github.com/mhibbins/heist
cd heist
python setup.py install
```

### From PyPI
```
pip install heist-hemiplasy
```

`heist` should be automatically added to your path.


## Usage
```
 _   _      ___ ____ _____
| | | | ___|_ _/ ___|_   _|
| |_| |/ _ \| |\___ \ | |
|  _  |  __/| | ___) || |
|_| |_|\___|___|____/ |_|
Hemiplasy Inference Simulation Tool
Version 0.3.0

Written by Mark Hibbins & Matt Gibson
Indiana University

usage: heist [-h] [-v] [-n] [-t] [-p] [-g] [-s] [-c] [-o] input

Tool for characterising hemiplasy given traits mapped onto a species tree

positional arguments:
  input                 Input NEXUS file

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Enable debugging messages to be displayed
  -n , --replicates     Number of replicates per batch
  -t , --threads        Number of threads for simulations
  -p , --mspath         Path to ms (if not in user path)
  -g , --seqgenpath     Path to seq-gen (if not in user path)
  -s , --mutationrate   Seq-gen mutation rate (default 0.05)
  -c , --CI             Optionally simulate at the upper ('upper') or lower
                        ('lower') bounds of the 95 % CI for the coalescent
                        conversion regression.
  -o , --outputdir      Output directory/prefix
```

## Input file

The input file is modified NEXUS format. A minimal example includes two trees (in newick format) and at least two dervived taxa set with the `set derived` command. If an outgroup is specified with `set outgroup`, the tree will be pruned to contain only taxa relevant to the simulation (i.e., the subclade containing derived taxa) + the outgroup. 

```
#NEXUS
begin trees;
tree tree_1 = (sp1:0.002,(sp2:0.001,((sp3:0.0004,sp4:0.0008)10.0:0.0005,(sp5:0.0006,sp6:0.0004)8.0:0.0004)15.0:0.0009)90.0:0.005);
tree tree_2 = (sp1:0.002,(sp2:0.001,((sp3:0.0004,sp4:0.0008)I1:0.0005,(sp5:0.0006,sp6:0.0004)I2:0.0004)I3:0.0009)I4:0.005)I5;
end;

begin hemiplasytool;
set derived taxon=sp2
set derived taxon=sp4
set derived taxon=sp6
end;
```

### Species trees
To run HeIST, you must supply two trees in Newick format. The first tree (must be named "tree_1") should have branch lengths in average substitutions per site and **branches must be labeled with concordance factors**. [IQTree](www.iqtree.org/doc/Concordance-Factor) can be used to do this. HeIST will use the concordance factors to convert to a tree in coalescent units (required for simulation).

The second tree (named "tree_2") should have the same branch lengths as tree_1, but have internal branch labels rather than concordance factors. These will allow the user to specify introgression events on internal branches.

> As of version 0.3.2, tree_2 must be supplied even if no introgression is being specified in the input.

> If your tree is already ultrametric and in coalescent units, you can supply this directly if you add the flag `set type coal` to the input file. In this case, tree_1 does not need to have concordance factors. 


### Traits

Set which taxa have the derived character by using 

```
set derived taxon="species in tree"
```


### Introgression

Introgression events can be defined by using

```
set introgression source="species in tree" dest="species in tree" prob=[float_value] timing=[float_value]
```
Note that timing must be specified in coalescent units. For this reason, we recommend first running your input tree through [`subs2coal`](#subs2coal)


### Conversion type parameter

Since ms requires input trees to be ultrametric, HeIST implements a tree smoothing step with two algorithm options: 

1. `ete3` -- redistributes branch lengths so that the distance from root to tip is the same, using the convert_to_ultrametric function in ete3 (default)
2. `extend` -- extends tip branches while preserving internal branch lengths

## Example

### NEXUS example

see `example/heist_example_input.txt`

```
#NEXUS
begin trees;
tree tree_1 = (sp1:0.002,(sp2:0.001,((sp3:0.0004,sp4:0.0008)10.0:0.0005,(sp5:0.0006,sp6:0.0004)8.0:0.0004)15.0:0.0009)90.0:0.005);
tree tree_2 = (sp1:0.002,(sp2:0.001,((sp3:0.0004,sp4:0.0008)I1:0.0005,(sp5:0.0006,sp6:0.0004)I2:0.0004)I3:0.0009)I4:0.005)I5;
end;

begin hemiplasytool;
set derived taxon=sp2
set derived taxon=sp4
set derived taxon=sp6
set introgression source=I1 taxon2=sp2 prob=0.05 timing=0.1
set conversion type=extend
end;
```

### HeIST call

see `example/heist_example.sh`

```
heist -v -n 100000 -t 16 -p ~/bin/ms -g ~/bin/seq-gen -o ./heist_example_output ./heist_example_input.txt
```



### Output

The above call will result in three output files being written: 

1. `heist_example_output.txt` contains the full summary
2. `heist_example_output.trees` contains observed gene trees in newick format
3. `heist_example_output_raw.txt` contains summary statistics in reduced format for merging multiple runs

```
### INPUT SUMMARY ###

Integer Code	Taxon Name
1:	sp1
2:	sp2
3:	sp3
4:	sp4
5:	sp5
6:	sp6
1:	I5
2:	I4
3:	I3
3:	I1
5:	I2

The species tree (smoothed, in coalescent units) is:
 (1:2.08154,(2:0.184423,((3:0.164423,4:0.164423)1:0.01,(5:0.164423,6:0.164423)1:0.01)1:0.01)1:1.89712);

Regression intercept: -0.24037164874509137
Regression slope: 424.7950852744482
X (newick internals): 0.0005,0.0004,0.0009,0.005
Y (coalescent internals): 0.01,0.01,0.01,1.8971199848858815
  ________________________________ 1
 |
_|                              ___ 2*
 |                             |
 |_____________________________| __ 3
                               ||
                               ||__ 4*
                               ||
                                |__ 5
                                |
                                |__ 6*

3 taxa have the derived state: 2, 4, 6

With homoplasy only, 3 mutations are required to explain this trait pattern (Fitch parsimony)

Introgression from taxon 3 into taxon 2 occurs at time 0.1 with probability 0.05

1.00e+05 simulations performed, using a mutation rate of 0.05

### RESULTS ###

127 loci matched the species character states

"True" hemiplasy (1 mutation) occurs 70 time(s)

Combinations of hemiplasy and homoplasy (1 < # mutations < 3) occur 52 time(s)

"True" homoplasy (>= 3 mutations) occurs 5 time(s)

127 loci have a discordant gene tree
0 loci are concordant with the species tree

6 loci originate from an introgressed history
121 loci originate from the species history

Distribution of mutation counts:

# Mutations	# Trees
On all trees:
1		70
2		52
3		5

On concordant trees:
# Mutations	# Trees

On discordant trees:
# Mutations	# Trees
1		70
2		52
3		5

Origins of mutations leading to observed character states for hemiplasy + homoplasy cases:

	Tip mutation	Internal branch mutation	Tip reversal
Taxa 2	0	52	0
Taxa 4	0	52	0
Taxa 6	0	52	0

### OBSERVED GENE TREES ###

                              _____ 4*
  ___________________________|
 |                           | ____ 2*
 |                           ||
_|                            |____ 6*
 |
 |          ______________________ 1
 |_________|
           |           ___________ 3
           |__________|
                      |___________ 5

  _________________________________ 1
 |
_|                               __ 3
 |                      ________|
 |                     |        |__ 5
 |_____________________|
                       |   _______ 6*
                       |__|
                          |      _ 2*
                          |_____|
                                |_ 4*

  ________________________________ 1
_|
 |   ______________________________ 3
 |__|
    |         _____________________ 5
    |________|
             |             ________ 6*
             |____________|
                          |      __ 2*
                          |_____|
                                |__ 4*

  _________________________________ 1
_|
 |                           ______ 5
 |__________________________|
                            | _____ 3
                            ||
                             |   __ 2*
                             |__|
                                |__ 4*
                                |
                                |__ 6*

                            _______ 6*
  _________________________|
 |                         | ______ 2*
 |                         ||
_|                          |______ 4*
 |
 |    _____________________________ 1
 |___|
     |                        _____ 3
     |_______________________|
                             |_____ 5

  _________________________________ 1
 |
_|                              ___ 3
 |                       ______|
 |                      |      |___ 5
 |______________________|
                        |     _____ 2*
                        |____|
                             | ___ 4*
                             ||
                              |___ 6*
...continued
```

From the above simulation, we can see that although the input species tree would have pointed to three independent origins of the arbitrary derived state (in taxa 2, 4, and 6), after we account for the possiblity of hemiplasy only 1 (or maybe 2) transitions are much more likely. 


## General guidelines for choosing the number of replicates 

Generally, the number of simulated loci with character states that match the observed distribution will be a small subset of the total number of loci. Therefore, it is typically necessary to simulate a large number of loci in order to observe a sufficient number of relevant cases. The precise number of loci to simulate will differ for each case, and will require some experimentation on the part of the user to come to an optimal value. We can provide some general guidelines to aid this exploration, however. Trees with fewer taxa and a higher specified mutation rate will require fewer simulations in order to observe relevant cases. The five-taxon test tree which uses a mutation rate value of 0.05, generates several hundred focal cases from 1x10^6 simulated loci. By contrast, the 15-taxon lizard phylogeny we analyze in our paper, which used a mutation rate of 0.001, required 1x10^10 simulations to observe a similar number of focal cases. Simulations of up to 1x10^7 loci are doable using the resources of a typical personal laptop, with memory use quickly becoming a limiting factor as the number of loci increases beyond this. We offer two approaches to aid with performance issues: 1) support for multiple processors, and 2) a module called “heistMerge” (see below) which combines the outputs from multiple independent runs. 


## Sub-modules

Once installed, two additional programs will be available at the command line: `newick2ms` and `subs2coal`.

### newick2ms

```
usage: newick2ms [-h] input

Tool for converting a newick string to ms-style splits. Note that this only
makes sense if the input tree is in coalescent units.

positional arguments:
  input       Input newick string file

optional arguments:
  -h, --help  show this help message and exit
```


### subs2coal

```
usage: subs2coal [-h] input

Tool for converting a newick string with branch lengths in subs/site to a
neewick string with branch lengths in coalescent units. Input requires gene or
site-concordancee factors as branch labels

positional arguments:
  input       Input newick string file

optional arguments:
  -h, --help  show this help message and exit
```

### heistMerge

```
usage: heistmerge [-h] [-d] [inputs [inputs ...]]

Merge output files from multiple HeiST runs. Useful for simulating large trees
by running multiple batch jobs.

positional arguments:
  inputs      Prefixes of output files to merge or a directory (supply -d flag
              as well)

optional arguments:
  -h, --help  show this help message and exit
  -d          Merge all files in a directory
```

`heistMerge` will write the merged output summary to standard out and create a new files `merged_trees.trees` which contains all observed focal gene trees.
