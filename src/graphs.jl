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

weight_matrix_to_vector(G::Graph{Directed},C) = map(e -> C[dst(e),src(e)], edges(G))
weight_matrix_to_vector(R::Field,G::Graph{Directed},C) = R.(weight_matrix_to_vector(G,C))

weight_vector_to_matrix(G::Graph{Directed},w) = weight_vector_to_matrix(MatElem,G,w)
weight_vector_to_matrix(R::Field,G::Graph{Directed},w) = weight_vector_to_matrix(MatElem,R,G,w)

weight_vector_to_matrix(::Type{Matrix},G::Graph{Directed},w) = Matrix(weight_vector_to_matrix(G,w))
weight_vector_to_matrix(::Type{Matrix},R::Field,G::Graph{Directed},w) = Matrix(weight_vector_to_matrix(R,G,w))

weight_vector_to_matrix(::Type{MatElem},G::Graph{Directed},w::Vector{<:Integer}) = weight_vector_to_matrix(MatElem,QQ,G,w) 
weight_vector_to_matrix(::Type{MatElem},G::Graph{Directed},w::Vector{<:FieldElem}) = weight_vector_to_matrix(MatElem,base_ring(w),G,w)
function weight_vector_to_matrix(::Type{MatElem},R::Field, G::Graph{Directed},w)
  length(w) != n_edges(G) && error("Weight vector does not correspond to edges of graph")
  C = identity_matrix(R,n_vertices(G))
  for (c,e) in zip(w,edges(G))
    C[dst(e),src(e)] = c
  end
  return C
end

indegree(G::Graph, v::Int) = inneighbors(G, v) |> length
outdegree(G::Graph, v::Int) = outneighbors(G, v) |> length


function transitively_closed_acyclic_graphs(n::Int)
    output = Graph{Directed}[]
    seen = Set{Int}()
    queue = Graph{Directed}[complete_dag(n)]
    while !is_empty(queue)
        G = popfirst!(queue)
        push!(output, G)
        E = transitive_reduction(G) |> edges |> collect
        for e in E
            rem_edge!(G, src(e), dst(e))
            graph_hash = Polymake.graph.canonical_hash(Oscar.pm_object(G))
            if graph_hash ∉ seen
              if edges(G) |> !is_empty
                H = graph_from_edges(Directed, edges(G), n)
                push!(queue, H)
              end
            end
            push!(seen, graph_hash)
            add_edge!(G, src(e), dst(e))
        end
    end
    push!(output, Graph{Directed}(n))
    return output
end

function vertex_to_edge_action(G::Graph, σ::PermGroupElem)
  if n_vertices(G) != degree(σ)
    throw("Permutation is not applicable to the graph")
  end

  E = edges(G) |> collect
  targets = map(E) do e
    σe = Edge(σ(src(e)), σ(dst(e)))
    findfirst(==(σe), E)
  end

  return perm(symmetric_group(n_edges(G)), targets)
end

