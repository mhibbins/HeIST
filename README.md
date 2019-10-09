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
python -m hemiplasytool -v -n 1000000 -p ~/bin/ms -g ~/bin/seq-gen -x 1 ./test/input_test_small.txt -o outtest.txt
```

### Output:
```
### INPUT SUMMARY ###

The species tree is (1:1.76011,(((3:0.990059,((8:0.33002,9:0.33002)1:0.33002,7:0.660039)1:0.33002)1:0.33002,(((10:0.440026,11:0.440026)1:0.440026,((12:0.33002,13:0.33002)1:0.33002,((15:0.220013,16:0.220013)1:0.220013,14:0.440026)1:0.220013)1:0.220013)1:0.220013,4:1.10007)1:0.220013)1:0.220013,(2:1.02673,(5:0.513364,6:0.513364)1:0.513364)1:0.513364)1:0.220013);

  ________________________________ 1
 |
 |              _________________ 3
 |             |
 |        _____|            _____ 8
 |       |     |      _____|
 |       |     |_____|     |_____ 9
 |       |           |
 |       |           |___________ 7
_|       |
 |       |                _______ 10
 |    ___|        _______|
 |   |   |       |       |_______ 11
 |   |   |       |
 |   |   |    ___|          _____ 12
 |   |   |   |   |    _____|
 |   |   |   |   |   |     |_____ 13
 |   |   |   |   |___|
 |   |   |   |       |        ___ 15
 |___|   |___|       |    ___|
     |       |       |___|   |___ 16
     |       |           |
     |       |           |_______ 14
     |       |
     |       |___________________ 4
     |
     |         __________________ 2
     |________|
              |          ________ 5
              |_________|
                        |________ 6

0 taxa have the derived state: 

The minimum number of mutations required to explain this trait pattern is 4


1000000 simulations performed
### OUTPUT SUMMARY ###

"True" hemiplasy (1 mutation) occurs 0 time(s)

Combinations of hemiplasy and homoplasy (1 < # mutations < 4) occur 0 time(s)

"True" homoplasy (>= 3 mutations) occurs 3 time(s)

3 loci have a discordant gene tree
0 loci are concordant with the species tree

0 loci originate from an introgressed history
3 loci originate from the species history

In cases with combinations of hemiplasy and homoplasy:

Taxon 2 reverted to the ancestral state 1 time(s).Taxon 4 reverted to the ancestral state 3 time(s).Taxon 5 reverted to the ancestral state 1 time(s).Taxon 14 reverted to the ancestral state 0 time(s).Taxon 15 reverted to the ancestral state 0 time(s).Taxon 16 reverted to the ancestral state 0 time(s).

### DETAILED OUTPUT ###

On concordant trees:
# Mutations	# Trees

On discordant trees:
# Mutations	# Trees
4		1
5		2

Derived mutation inheritance patterns for trees with fewer mutations than derived taxa:
	Term	Inherited from anc node
Taxa 2	1	2
Taxa 4	3	0
Taxa 5	1	2
Taxa 14	0	3
Taxa 15	0	3
Taxa 16	0	3

Of the replicates that follow species site pattern:
3 were discordant
0 were concordant

3 replicates matching the species pattern were from the species tree

### OBSERVED GENE TREES ###

  ________________________________ 1
 |
 |                 ______________ 10
 |           _____|
 |          |     |        ______ 12
 |          |     |_______|
 |       ___|             |______ 13
 |      |   |
 |      |   |  __________________ 11
_|      |   |_|
 |      |     |           _______ 14
 | _____|     |__________|
 ||     |                |    ___ 15
 ||     |                |___|
 ||     |                    |___ 16
 ||     |
 ||     |  ______________________ 3
 ||     |_|
 ||       |             _________ 9
  |       |____________|
  |                    |_________ 7
  |                    |
  |                    |_________ 8
  |
  |     __________________________ 4
  |____|
       |    _____________________ 6
       |___|
           |      _______________ 2
           |_____|
                 |_______________ 5

This topology occured 1 time(s)
  _______________________________ 1
 |
 |       ________________________ 11
_|      |
 |      |         _______________ 7
 |      |________|
 |      |        |         ______ 8
 |______|        |________|
        |                 |______ 9
        |
        |       _________________ 4
        |     _|
        |    | | ________________ 10
        |    | ||
        |    |  |     ___________ 12
        | ___|  |____|
        ||   |       |___________ 13
        ||   |
        ||   |____________________ 3
        ||   |
        ||   |         __________ 14
        ||   |________|
         |            |     _____ 15
         |            |____|
         |                 |_____ 16
         |
         |    ___________________ 5
         |___|
             |     ______________ 2
             |____|
                  |______________ 6

This topology occured 1 time(s)
         _________________________ 6
  ______|
 |      |        _________________ 2
 |      |_______|
 |              |_________________ 5
_|
 | _______________________________ 1
 ||
 ||                   ___________ 7
 ||     _____________|
  |    |             |___________ 9
  |    |
  |____|     _____________________ 3
       | ___|
       ||   |_____________________ 8
       ||
       ||    _____________________ 4
        |   |
        |   |         ____________ 12
        |___|  ______|
            | |      |____________ 13
            | |
            |_|         __________ 16
              |   _____|
              |  |     |  ________ 14
              |  |     |_|
              |__|       |________ 15
                 |
                 |  ______________ 10
                 |_|
                   |______________ 11

This topology occured 1 time(s)

...
```

See `outtest.txt` for full output.

![Mutation distribution](mutation_dist.png)


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