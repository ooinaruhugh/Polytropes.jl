using Oscar
import Oscar: RingElem 

@doc raw"""
    weighted_digraph_polyhedron(
        G::Graph{Directed}, 
        w::AbstractVector{<:RingElem};
        modulo_lineality=true
    )

The weighted digraph polyhedron $\mathrm{Q}(G)$ is defined by the linear inequalities 
$x_i - x_j \leq w_{ij}$ where $w_{ij}$ is the weight of edge $j\to i$.

The default behaviour is to return $\mathrm{Q}(G)$ modulo the subspace spanned by
$(1,1,\dots,1)$ which is controlled by `modulo_lineality`.
"""
function weighted_digraph_polyhedron(
        G::Graph{Directed}, 
        w::AbstractVector{<:RingElem};
        modulo_lineality=true
    )
    R = parent(w[1])
    A = sparse_matrix(R)
    if modulo_lineality
        for e in edges(G)
            row = [(dst(e)-1, 1)]
            if src(e) != 1
                push!(row, (src(e)-1, -1))
            end

            push!(A, sparse_row(R, row))
        end
    else
        for e in edges(G)
            push!(A, sparse_row(R, [(src(e), -1), (dst(e),1)]))
        end
    end

    return polyhedron(matrix(A), w)
end
