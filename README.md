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
  input                Input NEXUS file


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

The input file is modified NEXUS format. A minimal example includes a tree (in newick format) and at least two dervived taxa set with the `set derived` command.

```
#NEXUS
begin trees;
    tree tree_1 = (Lygosoma_sp_ACD7392:0.022135,(((Fojia_bumui_CCA15922:0.004971,((Papuascincus_stanleyanus_cf_sp2_CCA02483:0.002792,Papuascincus_sp_nov_CCA05980:0.003494)37.8:0.001417,Lipinia_pulchra_CCA03201:0.003904)41.9:0.002306)29.1:0.001868,(((Lipinia_albodorsale_CCA03500:0.003722,Lipinia_noctua_CCA16819:0.004171)20.2:0.001171,((Lobulia_elegans_CCA17485:0.001982,Lobulia_brongersmai_CCA03482:0.003384)31.5:0.001793,((Prasinohaema_sp_nov_CCA01005:0.002495,Prasinohaema_flavipes_CCA17083:0.002098)37:0.001049,Prasinohaema_prehensicauda_CCA17222:0.004007)33.3:0.001615)18.2:0.001013)4.37:1.87E-4,Prasinohaema_semoni_CCA00900:0.006197)6.49:3.78E-4)11.6:8.73E-4,(Prasinohaema_virens_CCA02661:0.008386,(Prasinohaema_sp_nov_CCA01623:0.008567,Lipinia_longiceps_CCA15866:0.006692)22.7:0.001413)17:0.001229)11.5:0.02697);
end;

begin hemiplasytool;
set derived taxon=Prasinohaema_sp_nov_CCA01005
set derived taxon=Prasinohaema_flavipes_CCA083
set derived taxon=Prasinohaema_prehensicauda_CCA222
set derived taxon=Prasinohaema_semoni_CCA00900
set derived taxon=Prasinohaema_virens_CCA02661
set derived taxon=Prasinohaema_sp_nov_CCA01623
end;

```

### Species tree

Species tree in newick format. Branch lengths must be in average substitutions per site and branches must be labeled with gene concordance factors.

### Traits

Set which taxa have the derived character by using 

```
set derived taxon="species in tree"
```


### Introgression

Introgression events can be defined by using

```
set introgression taxon1="species in tree" taxon2="species in tree" strength=[float_value]
```


## Example:
```
python -m hemiplasytool -v -n 1000000 -p ~/bin/ms -g ~/bin/seq-gen -x 50 ./test/input_test_small.txt -o outtest.txt
```

### Output:
```
### INPUT SUMMARY ###

Integer Code	Taxon Name
1:	Lygosoma_sp_ACD7392
2:	Prasinohaema_virens_CCA02661
3:	Prasinohaema_semoni_CCA00900
4:	Fojia_bumui_CCA15922
5:	Prasinohaema_sp_nov_CCA01623
6:	Lipinia_longiceps_CCA15866
7:	Lipinia_pulchra_CCA03201
8:	Lipinia_albodorsale_CCA03500
9:	Lipinia_noctua_CCA16819
10:	Papuascincus_stanleyanus_cf_sp2_CCA02483
11:	Papuascincus_sp_nov_CCA05980
12:	Prasinohaema_prehensicauda_CCA222
13:	Lobulia_elegans_CCA485
14:	Lobulia_brongersmai_CCA03482
15:	Prasinohaema_flavipes_CCA083
16:	Prasinohaema_sp_nov_CCA01005

The species tree (smoothed, in coalescent units) is:
 (((((((12:0.440047,(15:0.220024,16:0.220024)1:0.220024)1:0.220024,(13:0.330035,14:0.330035)1:0.330035)1:0.220024,(8:0.440047,9:0.440047)1:0.440047)1:0.220024,3:1.10012)1:0.220024,((7:0.660071,(10:0.330035,11:0.330035)1:0.330035)1:0.330035,4:0.990106)1:0.330035)1:0.220024,((5:0.513388,6:0.513388)1:0.513388,2:1.02678)1:0.513388)1:0.220024,1:1.76019);

                          ______ 12*
                      ___|
                     |   |    ___ 15*
                     |   |___|
                  ___|       |___ 16*
                 |   |
                 |   |     _____ 13
              ___|   |____|
             |   |        |_____ 14
             |   |
          ___|   |       _______ 8
         |   |   |______|
         |   |          |_______ 9
         |   |
      ___|   |___________________ 3*
     |   |
     |   |           ___________ 7
     |   |     _____|
     |   |    |     |      _____ 10
     |   |____|     |_____|
  ___|        |           |_____ 11
 |   |        |
 |   |        |_________________ 4
 |   |
 |   |                  ________ 5*
_|   |         ________|
 |   |________|        |________ 6
 |            |
 |            |__________________ 2*
 |
 |______________________________ 1

6 taxa have the derived state: 12, 15, 16, 3, 5, 2

The minimum number of mutations required to explain this trait pattern is 4


100000000 simulations performed
### OUTPUT SUMMARY ###

"True" hemiplasy (1 mutation) occurs 11 time(s)

Combinations of hemiplasy and homoplasy (1 < # mutations < 4) occur 42 time(s)

"True" homoplasy (>= 3 mutations) occurs 24 time(s)

77 loci have a discordant gene tree
0 loci are concordant with the species tree

0 loci originate from an introgressed history
77 loci originate from the species history

In cases with combinations of hemiplasy and homoplasy:

Taxon 2 mutated to the derived state 10 time(s), and inherited it from an ancestral population 53 time(s)
Taxon 3 mutated to the derived state 6 time(s), and inherited it from an ancestral population 57 time(s)
Taxon 5 mutated to the derived state 13 time(s), and inherited it from an ancestral population 50 time(s)
Taxon 12 mutated to the derived state 2 time(s), and inherited it from an ancestral population 61 time(s)
Taxon 15 mutated to the derived state 1 time(s), and inherited it from an ancestral population 62 time(s)
Taxon 16 mutated to the derived state 0 time(s), and inherited it from an ancestral population 63 time(s)
Taxon 1 reverted to the ancestral state 0 time(s).Taxon 4 reverted to the ancestral state 0 time(s).Taxon 6 reverted to the ancestral state 1 time(s).Taxon 7 reverted to the ancestral state 0 time(s).Taxon 8 reverted to the ancestral state 0 time(s).Taxon 9 reverted to the ancestral state 0 time(s).Taxon 10 reverted to the ancestral state 0 time(s).Taxon 11 reverted to the ancestral state 0 time(s).Taxon 13 reverted to the ancestral state 0 time(s).Taxon 14 reverted to the ancestral state 0 time(s).

### DETAILED OUTPUT ###

On concordant trees:
# Mutations	# Trees

On discordant trees:
# Mutations	# Trees
1		11
2		21
3		21
4		16
5		5
6		2
7		1

Derived mutation inheritance patterns for trees with fewer mutations than derived taxa:
	Term	Inherited from anc node
Taxa 2	10	53
Taxa 3	6	57
Taxa 5	13	50
Taxa 12	2	61
Taxa 15	1	62
Taxa 16	0	63
Taxa 1	0	1
Taxa 4	0	1
Taxa 6	1	0
Taxa 7	0	1
Taxa 8	0	1
Taxa 9	0	1
Taxa 10	0	1
Taxa 11	0	1
Taxa 13	0	1
Taxa 14	0	1

Of the replicates that follow species site pattern:
77 were discordant
0 were concordant

77 replicates matching the species pattern were from the species tree

### OBSERVED GENE TREES ###

                  _______________ 3
             ____|
            |    |   _____________ 12
            |    |__|
            |       |       ______ 15
  __________|       |______|
 |          |              |______ 16
 |          |
 |          |    _________________ 2
 |          |___|
_|              |_________________ 5
 |
 |     ___________________________ 1
 |    |
 |    |           ________________ 11
 |____| _________|
      ||         |  _____________ 4
      ||         |_|
      ||           |_____________ 7
      ||           |
       |           |_____________ 10
       |
       |    ______________________ 6
       |   |
       |___|               ______ 8
           |    __________|
           |   |          |______ 9
           |___|
               |           _______ 13
               |__________|
                          |_______ 14

This topology occured 1 time(s)
                        _________ 3
                    ___|
                   |   |  _______ 15
                   |   |_|
                   |     |    ____ 12
  _________________|     |___|
 |                 |         |____ 16
 |                 |
 |                 |   __________ 6
 |                 |__|
 |                    |  ________ 2
_|                    |_|
 |                      |________ 5
 |
 |                _______________ 1
 |               |
 |               |       ________ 4
 |_______________|   ___|
                 |  |   |________ 7
                 |  |
                 |__|       _____ 10
                    |______|
                    |      |_____ 11
                    |
                    |     _______ 13
                    |____|
                         |_______ 14
                         |
                         |    ___ 8
                         |___|
                             |___ 9

This topology occured 1 time(s)
                         ________ 2
                    ____|
                   |    |________ 5
  _________________|
 |                 |    __________ 3
 |                 |___|
 |                     |    _____ 12
 |                     |___|
_|                         |   __ 15
 |                         |__|
 |                            |__ 16
 |
 |               ________________ 6
 |              |
 |______________|  ______________ 1
                | |
                |_|   ___________ 7
                  |  |
                  |__|      _____ 10
                     | ____|
                     ||    |_____ 11
                     ||
                      |__________ 4
                      |
                      |        __ 13
                      |   ____|
                      |  |    |__ 14
                      |__|
                         |    ___ 8
                         |___|
                             |___ 9


...
```

See `lizard_tree_output.txt` for full output.

![Mutation distribution](mutation_dist_lizard.png)


## Sub-modules

Once installed, two additional programs will be available at the command line: `newick2ms` and `subs2coal`.

### `newick2ms`

```
usage: newick2ms [-h] input

Tool for converting a newick string to ms-style splits. Note that this only
makes sense if the input tree is in coalescent units.

positional arguments:
  input       Input newick string file

optional arguments:
  -h, --help  show this help message and exit
```


### `subs2coal`

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