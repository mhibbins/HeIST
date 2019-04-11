#/usr/bin/bash

### Test run for a 6-taxon ultrametric tree (figure out tip branches later) 

#call ms 
msdir/ms 6 300000 -T -I 6 1 1 1 1 1 1 -ej 6 2 1 -ej 3 3 2 -ej 1.5 5 3 -ej 1.25 6 5 -ej 1 4 3 | tail -n +4 | grep -v // >test_sptree_ultra.txt 

#call seq-gen 
Seq-Gen/source/seq-gen -m HKY -l 1 -s 0.01 -wa <test_sptree_ultra.txt >test_sptree_ultra_seqs.txt

