#!/bin/bash

N=$1

count=`ls -1 data/combinatorial-types/$N/*.topcom | wc -l`

echo

for f in data/combinatorial-types/$N/*.topcom; do
  mpirun -np $(nproc) mptopcom --regular <$f 2>/dev/null \
    | pv -cl --name $(basename $f) | xz > data/combinatorial-types/$N/$(basename $f .topcom).regular.out.xz
  printf "{$f}\n" 
done | pv -cl --size $count --name Total > /dev/null

