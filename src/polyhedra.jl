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
        w::Union{AbstractMatrix{T},MatElem{T},AbstractVector{T}};
        modulo_lineality=true
) where {T<:RingElem}
    #R = base_ring(w)
    R = parent(first(w))
    A = if modulo_lineality
        fundamental_polytope(Matrix, G, R)[1:end-1, 2:end]
    else 
        fundamental_polytope(Matrix, G, R)[1:end-1, :]
    end

    return polyhedron(A, collect(w) |> vec)
end

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
        w::Union{AbstractMatrix{T},AbstractVector{T}};
        modulo_lineality=true
) where {T<:Real}
    A = if modulo_lineality
        fundamental_polytope(Matrix, G, T)[1:end-1, 2:end]
    else 
        fundamental_polytope(Matrix, G, T)[1:end-1, :]
    end

    return polyhedron(A, collect(w) |> vec)
end

function tropical_ball(
    center::Union{AbstractMatrix{T},AbstractVector{T}}, 
    radius::T; 
    modulo_lineality=true
) where {T<:Real}
    G = complete_directed_graph(length(center))
    A = fundamental_polytope(Matrix, G, T)[1:end-1, :]
    b = radius .+ A * center

    return weighted_digraph_polyhedron(G, b; modulo_lineality=modulo_lineality)
end
