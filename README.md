# HemiplasyTool

## Authors:
Matt Gibson (gibsomat@indiana.edu)  
Mark Hibbins (mhibbins@indiana.edu)

## Dependencies:
* ms  
* seq-gen  
* biopython
* scipy
* numpy


## Usage
```
usage: hemiplasytool.py [-h] [-v] [-n] [-x] [-p] [-g] [-o] splits traits tree

Calculate the probability that convergent trait patterns are due to hemiplasy

positional arguments:
  splits              Split times file, ordered from oldest to newest. In
                      units of 4N generations.
  traits              Traits file
  tree                Species topology in Newick format on one line.

optional arguments:
  -h, --help          show this help message and exit
  -v, --verbose       Enable debugging messages to be displayed
  -n , --replicates   Number of replicates per batch
  -x , --batches      Number of batches
  -p , --mspath       Path to ms
  -g , --seqgenpath   Path to seq-gen
  -o , --outputdir    Output directory
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
python hemiplasytool/hemiplasytool.py -v -n 1000000 -p ~/bin/ms -g ~/bin/seq-gen -x 2 test/test_splits.txt test/test_traits.txt test/test_topology.txt
```

### Output:
```
Of the replicates that follow species site pattern:
324 were discordant
229 were concordant

Of the replicates that did not follow the species site pattern:
1137961 were discordant
861486 were concordant

Odds ratio: 1.0711008741372166
P-val: 0.4397051096565977

On concorant trees:
# Mutations	# Trees
1		5
2		35
3		175
4		13
5		1

On discordant trees:
# Mutations	# Trees
1		9
2		43
3		250
4		17
5		4
6		1

DEBUG:root:Plotting...

Time elapsed: 324.11109495162964 seconds
```

![Mutation distribution](mutation_dist.png)