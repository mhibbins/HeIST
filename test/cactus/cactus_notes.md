# Notes for application of hemiplasytool to cactus phylogeny from Copetti et al. 2017

Here we present some notes relating to our application of hemiplasytool
to the tree presented in Figure 2 of 

Copetti D., Burquez A., Bustamante E., ... , Sanderson M.J. 2017. 
Extensive gene tree discordance and hemiplasy shaped the genomes of 
North American columnar cacti. *PNAS* 114(45): 12003 - 12008. 

## Taxa labels 

We used the following labelling conventions for the tree:

1: *Pereskia humboldtii*
2: *Stenocereus thurberi*
3: *Lophocereus schottii*
4: *Carnegiea gigantea*
5: *Pacycereus pringlei*

## Demographic parameters

We used the following demographic parameter estimates, adapted from 
the paper, to estimate split times in the required units of 4N 
generations: 

Generation time: 45 years
Effective population size: 30,000

1/2 split: (27 mya / 45) / (4 * 30,000) = 5 
2/3 split: (8 mya / 45) / (4 * 30,000) = 1.5
3/4 split: (5 mya / 45) / (4 * 30,000) = 0.9
4/5 split: (3 mya / 45) / (4 * 30,000) = 0.5

We also specify the mutation rate in units of per site per 4N 
generations: 

6.0 x 10^-8 per site per generation * 4 * 30,000 = 7.2 x 10^-3

## Traits 

The authors study patterns of amino acid sequence variation, but do
not report any particular re-occuring patterns. Therefore, we simulate
several different configurations of derived states that could plausibly
arise due to hemiplasy. 

## Introgression

The authors infer two introgression events from their network analysis.
However, one of these (the green reticulation in Figure 2) involves an
unsampled ghost lineage, and so we have ignored it. The blue reticulation
is specified from 2 into 5, at a time of 0.15 with a proportion of 6.5%. 
 
