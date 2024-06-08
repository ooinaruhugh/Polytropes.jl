using Oscar
import Oscar: ZZRingElem, RayVector, QQFieldElem

function polytrope(G::Graph{Directed}, w::AbstractVector{<:RingElem})
  R = parent(w[1])
  n = number_of_vertices(G)
  A = sparse_matrix(R)

  for e in edges(G)
    push!(A, sparse_row(R, [(src(e), -1), (dst(e),1)]))
  end

  return polyhedron(matrix(A), w)
end

@doc raw"""
    weighted_digraph_polyhedron(G::Graph{Directed}, w::AbstractVector{<:RingElem})

The weighted digraph polyhedron $Q(G)$ is defined by the linear inequalities 
$x_i - x_j \leq w_{ij}$ where $w_{ij}$ is the weight of edge $j\to i$.
"""
function weighted_digraph_polyhedron(
        G::Graph{Directed}, 
        w::AbstractVector{<:RingElem}
    )
    # This is already happening in the tropical projective torus
    R = parent(w[1])

    n = number_of_vertices(G) - 1
    A = sparse_matrix(R)

    for e in edges(G)
        row = [(dst(e)-1, 1)]
        if src(e) != 1
            push!(row, (src(e)-1, -1))
        end

        push!(A, sparse_row(R, row))
    end

    return polyhedron(matrix(A), w)
end
