using ProgressBars

@assert length(ARGS) == 1 "Usage: calculate-stratum.jl MRDI_FILE"
input_file, = ARGS

include("compute-stratum.jl")

compute_stratum(input_file)

