# Polytropes.jl

This repository contains a Julia module for computations on *polytropes* using OSCAR.


## Enumerating combinatorial types of polytropes
First, generate the point configurations corresponding to each DAG on n nodes. 
Run the following command from the top-level of the source tree to generate the files for $n=4$.
```sh
julia --project examples/with-mptopcom.jl 4
```
Then, run `mptopcom` for each file in `tmp/<n>/`. The following command automates this and calculates the total.
```sh
export N=4; for file in tmp/$N/*.topcom; do mptopcom1 --central < $file 2> /dev/null | wc -l; done | awk '{s+=$1} END {print s}'
```
