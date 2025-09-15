#!/bin/bash

N=$1

time (for file in tmp/$N/*.topcom; do mptopcom1 --central --regular < $file 2> /dev/null | wc -l ; printf "{$file}\n" 1>&2; done) > $N.counts
awk '{s+=$1} END {print s+1}' $N.counts
