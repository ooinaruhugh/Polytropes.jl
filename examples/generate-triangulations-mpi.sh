#!/bin/bash

folder=$1

count=`ls -1 $folder/*.topcom | wc -l`

echo

for f in $folder/*.topcom; do
  mpirun -np $(nproc) mptopcom --central --regular <$f 2>/dev/null \
    | pv -cl --name $(basename $f) | xz > $folder/$(basename $f .topcom).out.xz
  printf "{$f}\n" 
done | pv -cl --size $count --name Total > /dev/null

