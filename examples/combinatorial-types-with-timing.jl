using Polytropes

counts = Int[]
time_elapsed = Float64[]

for i in 1:6
    t = @elapsed begin
        n = i |> complete_dag |> enumerate_satisfying_assignments |> length
    end

    push!(counts, n)
    push!(time_elapsed, t)
end

n_t = zip(counts, time_elapsed) |> collect