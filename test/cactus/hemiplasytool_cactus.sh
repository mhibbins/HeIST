#!/bin/bash

cd /N/dc2/scratch/mhibbins/hemiplasytool/
module unload python/2.7.13
module load python/3.6.1

for i in {1..5}
do
python -m hemiplasytool -v -n 500000 -p msdir/ms -g Seq-Gen/source/seq-gen -x 1 -s 0.0072 -o test/cactus/  test/cactus/cactus_input_$i.txt > test/cactus/cactus_output_$i.txt
done

