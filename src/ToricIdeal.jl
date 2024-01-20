using Oscar
using Bijections

import Base: Vector
import Oscar: toric_ideal, Graph, Undirected, Directed
# import Graphs: vertices, Edge

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
    E = sort(edges(G) |> collect; by=x->x.target)

    return Bijection(Dict(zip(E, 1:length(E))))
end

function variables(G::Graph{Directed})
    E = sort(edges(G) |> collect, by=x->x.target)
    return map(e->"e$(e.source)$(e.target)", E)
end

@doc raw"""
    edge_ring((G::Graph{Directed})

Returns for a given graph `G` the polynomial ring in variables $e_ij$ for every edge $i\to j$.

# Examples
```jldoctest
julia> complete_dag(3) |> edge_ring
Multivariate polynomial ring in 3 variables e12, e13, e23
  over rational field

julia> complete_dag(4) |> edge_ring
Multivariate polynomial ring in 6 variables e12, e13, e23, e14..., e34
  over rational field
```
"""
function edge_ring(G::Graph{Directed})
    #TODO: Change to return also vars
    return first(polynomial_ring(QQ, variables(G)))
end


function tropical_edge_ring(G::Graph{Directed}; minmax=min)
    #TODO: Change to return also vars
    return first(polynomial_ring(tropical_semiring(minmax), variables(G)))
end

@doc raw"""
    design_matrix((G::Graph{Directed})

Returns for a given graph `G` the design matrix, that is, the constraint matrix 
associated to the all-pairs shortest-paths problem on $G$.

# Examples
```jldoctest
julia> complete_dag(3) |> design_matrix
[ 1    1    0]
[-1    0    1]
[ 0   -1   -1]

julia> complete_dag(4) |> design_matrix
[ 1    1    0    1    0    0]
[-1    0    1    0    1    0]
[ 0   -1   -1    0    0    1]
[ 0    0    0   -1   -1   -1]
```
"""
function design_matrix(G::Graph{Directed})
    #TODO: I think this is called a presentation matrix (cf. semigroup presentations)
    I = indices(G)
    A = zero_matrix(ZZ, nv(G), length(I))

    for j in 1:nv(G)
        for i in inneighbors(G,j)
            A[j, get(I, Edge(i,j), 0)] = -1
        end
        for i in outneighbors(G,j)
            A[j, get(I, Edge(j,i), 0)] = 1
        end
    end

    return A
end

@doc raw"""
    toric_ideal((G::Graph{Directed})

Returns for a given graph `G` the toric ideal of the design matrix in the edge ring.

# Examples
```jldoctest
julia> complete_dag(3) |> toric_ideal
ideal(e12*e23 - e13)

julia> complete_dag(4) |> toric_ideal
ideal(e12*e24 - e14, e13*e34 - e14, e23*e34 - e24)
```
"""
function toric_ideal(G::Graph{Directed})
    R = edge_ring(G)
    MB = design_matrix(G) |> markov_basis

    return binomial_exponents_to_ideal(R, MB)
end

# function parse_gfan_output(s::String)
#     unnecessary_stuff = ['{', '}', '\n', ',', ' ']

#     cleaned_bases = s |>
#         (x -> split(x, "{{")[2]) |> 
#         (x -> strip(x, unnecessary_stuff)) |> 
#         (x -> split(x, "}\n")) .|>
#         (x -> strip(x, unnecessary_stuff)) .|>
#         (x -> split(x, ","))
        
#     output = map(x -> map(y -> strip(y), x), cleaned_bases)
#     initial_terms = map(y -> map(x -> split(x, ['+', '-'])|>first, y), output)

#     return output, initial_terms
# end

# function construct_poly(B::Vector{String})
#     return map(f -> eval(Meta.parse(f)), B)
# end

# function construct_poly(B::Vector{SubString{String}})
#     return construct_poly(String.(B))
# end

# function gfan(I::Ideal)
#     basis = gens(I)
#     vars = gens(base_ring(I))

#     input = "Q[" * join(vars, ",") * "]{" * join(basis, ",") * "}"

#     # println(input)

#     output = read(pipeline(`gfan_bases`, stdin=IOBuffer(input)), String)
#     parsed, parsed_initials = parse_gfan_output(output)
#     # bases = construct_poly.(parsed)
#     # initials = construct_poly.(parsed_initials)

#     return parsed, parsed_initials
#     # return bases, initials
# end