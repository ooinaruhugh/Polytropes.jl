#!/bin/bash

N=$1

time (for file in tmp/$N/*.topcom; do mpirun -np $(nproc) mptopcom --central --regular < $file 2> /dev/null | wc -l ; printf "{$file}\n" 1>&2; done) > $N.counts
awk '{s+=$1} END {print s+1}' $N.counts
