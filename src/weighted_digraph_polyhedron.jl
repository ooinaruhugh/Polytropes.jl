using Oscar

Vector{ZZRingElem}(v::RayVector{QQFieldElem}) = ZZ.(convert(Vector, v))

function weighted_digraph_polyhedron(G::Graph{Directed}, w)
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