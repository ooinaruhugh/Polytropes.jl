using Oscar
import Oscar: ZZRingElem, RayVector, QQFieldElem

Vector{ZZRingElem}(v::RayVector{QQFieldElem}) = ZZ.(convert(Vector, v))

function polytrope(G::Graph{Directed}, w::AbstractVector)
  n = number_of_vertices(G)
  A = sparse_matrix(QQ)

  for e in edges(G)
    push!(A, sparse_row(QQ, [(src(e), -1), (dst(e),1)]))
  end

  return polyhedron(MatElem(A), w)
end

function weighted_digraph_polyhedron(G::Graph{Directed}, w)
    # This is already happening in the tropical projective torus
    n = number_of_vertices(G) - 1
    A = []

    for e in edges(G)
        row = zero(1:n)
        if src(e) != 1
            row[src(e)-1] = -1
        end
        # if dst(e) != 1
            row[dst(e)-1] = 1
        # end

        push!(A, row)
    end

    A = reduce(hcat, A)'

    return polyhedron(A, w)
end
