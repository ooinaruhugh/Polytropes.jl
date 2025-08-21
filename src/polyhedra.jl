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
        fundamental_polytope(Matrix, G, T)[2:end, 2:end]
    else 
        fundamental_polytope(Matrix, G, T)[2:end, :]
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
    
    return (A[1:m+1,n+1:end] - A[1:m+1,1:n])[[end, 1:m...],:]
end
fundamental_polytope(G::Graph, R=ZZ) = fundamental_polytope(Matrix, G, R) |> convex_hull

function subdivision_of_fundamental_polytope(
    G::Graph{Directed}, w::AbstractVector{T}; project=true, acyclic=true
) where {T}
  A = fundamental_polytope(Matrix, G)

  if project
    S = subdivision_of_points(A[:, 2:end], [0,w...])
    if acyclic
      I = reduce(vcat,Oscar.pm_object(S).POLYHEDRAL_COMPLEX.MAXIMAL_POLYTOPES_INCIDENCES)
      subI = I[filter(x -> !I[x,1], 1:n_rows(I)), 2:end]
      return subdivision_of_points(A[2:end,2:end-1], subI)
    else
      return S
    end
  else
    return subdivision_of_points(A,[w...,0])
  end
end

function secondary_fan(A)
  #P = Polymake.polytope.PointConfiguration(POINTS=homogenize(A))
  P = Polymake.polytope.PointConfiguration(POINTS=A)
  return Polymake.fan.secondary_fan(P) |> polyhedral_fan
end
