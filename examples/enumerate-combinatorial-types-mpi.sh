#!/bin/bash

N=$1

count=`ls -1 $N/*.topcom | wc -l`

echo

for f in $N/*.topcom; do
  mpirun -np $(nproc) mptopcom --regular <$f 2>/dev/null \
    | pv -cl --name $(basename $f) | xz > $N/$(basename $f .topcom).regular.out.xz
  printf "{$f}\n" 
done | pv -cl --size $count --name Total > /dev/null

