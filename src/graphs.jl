using Oscar
using LinearAlgebra
import Oscar: Graph, Undirected, Directed

@doc raw"""
    complete_dag(n::Int64)

Returns the complete acyclic directed graph on $n$ nodes.

# Examples
```jldoctest
julia> complete_dag(3)
Directed graph with 3 nodes and the following edges:
(1, 2)(1, 3)(2, 3)
```
"""
function complete_dag(n)
    G = Graph{Directed}(n);

    for j in 1:n
        for i in (j+1):n
            add_edge!(G, j, i)
        end
    end

    return G
end

function complete_directed_graph(n)
    return graph_from_adjacency_matrix(Directed, ones(Bool, n, n) - I(n))
end

@doc raw"""
    indices((G::Graph{Directed})

Returns the dictionary mapping edges to the index of the corresponding variable 
in the edge ring and the design matrix.

# Examples
```jldoctest
julia> complete_dag(3) |> indices
Bijection Dict{Edge, Int64} with 3 entries:
  Edge(1, 2) => 1
  Edge(1, 3) => 2
  Edge(2, 3) => 3
```
"""
function indices(G::Graph{Directed})
    E = edges(G)

    return Dict(zip(E, 1:length(E)))
end

function opposite_graph(K::Graph{Directed})
    return graph_from_edges(Directed, [[dst(e), src(e)] for e in edges(K)], n_vertices(K))
end

function essential_edges(G::Graph{Directed})
    Gt = transitive_reduction(G)

    return GC.@preserve Gt (Gt |> edges |> collect)
end


indegree(G::Graph, v::Int) = inneighbors(G, v) |> length
outdegree(G::Graph, v::Int) = outneighbors(G, v) |> length

function root_polytope(::Type{Matrix}, G::Graph, R=ZZ)
    n = n_vertices(G)
    s = edges(G) .|> src
    t = edges(G) .|> dst
    
    return R.(hcat(I[[s...,1:n...],1:n], I[[t...,1:n...],1:n]))
end
root_polytope(G::Graph, R=ZZ) = root_polytope(Matrix, G, R) |> convex_hull

function fundamental_polytope(::Type{Matrix}, G::Graph, R=ZZ)
    A = root_polytope(Matrix, G, R)
    
    n = n_vertices(G)
    m = n_edges(G)
    
    return A[1:m+1,n+1:end] - A[1:m+1,1:n]
end
fundamental_polytope(G::Graph, R=ZZ) = fundamental_polytope(Matrix, G, R) |> convex_hull