# Polytropes.jl

This repository contains a Julia module for computations on *polytropes* using OSCAR.


## Enumerating combinatorial types of polytropes
This can be done by calling the following script from the top-level directory of this repo:
```sh
examples/enumerate-with-mptopcom.sh <n>
```

This generates the necessary point configurations corresponding to each DAG on $n$ nodes, and then feeds them to `mptopcom`.

One can also generate the point configurations separately by calling the following command.
```sh
julia --project examples/with-mptopcom.jl 4
```
Then, run `mptopcom` for each file in `tmp/<n>/`. The following command automates this and calculates the total.
```sh
export N=4; for file in tmp/$N/*.topcom; do mptopcom1 --central < $file 2> /dev/null | wc -l; done | awk '{s+=$1} END {print s}'
```
