using Oscar
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
    G = Graph{Directed}(n);

    for j in 1:n
        for i in 1:n
            i != j && add_edge!(G, j, i)
        end
    end

    return G
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
    #E = edges_by_target(G)
    E = edges(G)

    return Bijection(Dict(zip(E, 1:length(E))))
end

function variable_labels(G::Graph{Directed})
    #E = edges_by_target(G)
    E = edges(G)
    return map(e->"e$(e.source)$(e.target)", E)
end

function edges_by_target(G::Graph{Directed})
    return sort(edges(G) |> collect; by=x->x.target)
end

function edge_of_gen(G::Graph{Directed}, x::QQMPolyRingElem)
    return edges_by_target(G)[findfirst(y->y==x, edge_ring((G)) |> gens)]
end

function opposite_graph(K::Graph{Directed})
    return graph_from_edges(Directed, [[dst(e), src(e)] for e in edges(K)], n_vertices(K))
end

indegree(G::Graph, v::Int) = inneighbors(G, v) |> length
outdegree(G::Graph, v::Int) = outneighbors(G, v) |> length

function vertices_of_newton_polytope(G::Graph{Directed})
    C = Vector{PointVector}[]
    n = n_vertices(G)

    for v in 1:n
        Cv = [
            point_vector([i==u for i in  1:n])
            for u in [v,outneighbors(G, v)...]
        ]

        push!(C,Cv)
    end

    return C
end

function cayley_embedding_of_dual_vertices(G::Graph{Directed})
    C = vertices_of_newton_polytope(G)
    return cayley_embedding(C...)
end
