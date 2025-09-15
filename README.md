# Polytropes.jl

This repository contains a Julia module for computations on *polytropes* using OSCAR.

## Enumerating triangulations for fundamental polytopes
This can be done by calling the following scripts from the top-level directory of this repo:
```sh
julia --project examples/with-mptopcom.jl <n>
examples/enumerate-with-mptopcom.sh <n>
```

This generates the necessary point configurations corresponding to each DAG on $n$ nodes, and then feeds them to `mptopcom`.

## Enumerating combinatorial types of polytropes
This is done by calculating from above triangulations the corresponding secondary cones and counting the number
of cones (also the lower-dimensional ones) in the resulting fans.

For this, generate the triangulations with `examples/with-mptopcom.jl` as above. Then call the following scripts:
```sh
examples/generate-triangulations-mpi.sh tmp/<n>
julia examples/parse-triangulations.batch.jl tmp/<n>
julia examples/compute-stratum.batch.jl tmp/<n>
```
