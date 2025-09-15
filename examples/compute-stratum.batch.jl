using ProgressBars

@assert length(ARGS) == 1 "Usage: calculate-stratum.jl folder"
folder, = ARGS
files = readdir(folder; join=true) |> filter(contains(r"\d.mrdi"))

for input_file in ProgressBar(files)
  compute_stratum(input_file)
end
