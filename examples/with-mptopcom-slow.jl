using Oscar
using Polymake
using Polytropes
using Combinatorics

using ProgressBars

if length(ARGS) < 1
  n = 4
else
  n = parse(Int,ARGS[1])
end

if !isdir("tmp/$n")
  mkpath("tmp/$n")
end

E = complete_dag(n) |> edges |> collect

seen = Set(Int[])
GG = Graph{Directed}[]
for E in powerset(E) |> ProgressBar
  if !isempty(E) 
    local G = graph_from_edges(Directed, E, n)
    h = Polymake.graph.canonical_hash(G)

    if !(h in seen)
      push!(GG, G)
      push!(seen, h)
    end
  end
end

for (i,G) in ProgressBar(enumerate(GG))
  P = Polytropes.fundamental_polytope(G)
  pP = P |> Oscar.pm_object |> Polymake.polytope.project_full
  V = pP.VERTICES
  aut = pP.VERTICES_IN_FACETS |> Polymake.graph.automorphisms .|> last |> filter(==(0)∘first)

  open("tmp/$n/$i.topcom", "w+") do io
    println(io, "[")
      for i in 1:size(V, 1)
        println(io, Vector{Int}(V[i,:]), ",")
      end
    println(io, "]")

    if !isempty(aut)
      println(io, "[")
        for a in aut
          println(io, Vector{Int}(a), ",")
        end
      println(io, "]")
    end
  end
end

