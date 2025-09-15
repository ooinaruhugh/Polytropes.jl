#!/bin/bash

N=$1

#if [ ! -d tmp/$N ]; then
#  mkdir -p tmp/$N
#  julia --project examples/with-mptopcom.jl $N
#fi

count=`ls -1 data/combinatorial-types/$N/*.topcom | wc -l`

echo

for f in data/combinatorial-types/$N/*.topcom; do
  mpirun -np $(nproc) mptopcom --central --regular <$f 2>/dev/null \
    | pv -cl --name $(basename $f) | xz > data/combinatorial-types/$N/$(basename $f .topcom).out.xz
  printf "{$f}\n" 
done | pv -cl --size $count --name Total > /dev/null

