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

G = complete_dag(n)
E = edges(G) |> collect
PE = powerset(E) |> collect

GG = unique(last, map(PE |> filter(!isempty) |> ProgressBar) do E
   G = graph_from_edges(Directed, E, n)
   G, Polymake.graph.canonical_hash(G)
 end |> ProgressBar) .|> first

P = ProgressBar(Polytropes.fundamental_polytope.(GG)) .|> Oscar.pm_object .|> Polymake.polytope.project_full
#A = map(P) do p
A = map(ProgressBar(P)) do p
  V = p.VERTICES
  aut = p.VERTICES_IN_FACETS |> Polymake.graph.automorphisms .|> last |> filter(==(0)∘first)

  Matrix{Int}(V), Vector{Int}.(aut)
end

counts = Int[]

#for (i,(V,aut)) in enumerate(A)
for (i,(V,aut)) in ProgressBar(enumerate(A))
  open("tmp/$n/$i.topcom", "w+") do io
    println(io, "[")
      for i in 1:size(V, 1)
        println(io, V[i,:], ",")
      end
    println(io, "]")

    if !isempty(aut)
      println(io, "[")
        for a in aut
          println(io, a, ",")
        end
      println(io, "]")
    end
  end
end
#GG = transitively_closed_acyclic_graphs(n)
