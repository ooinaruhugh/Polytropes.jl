using Oscar
import Oscar: minkowski_sum

@doc raw"""
    minkowski_sum(A::AbstractVector{<:PointVector})

Return the Minkowski sum $\{ \sum_{i=1}^n a_i \mid a_i\in A_i \}$ of
the point sets in `A`.
"""
function minkowski_sum(A::AbstractVector{<:PointVector}...)
    return Iterators.product(A...) .|> sum |> unique
end

function points_to_total_indices(A::AbstractVector{<:PointVector}...)
    l = 0
    combined_indices = Vector{Int}[]
    for (j,a) in enumerate(A)
        push!(combined_indices, (enumerate(a).|>first).+l)
        l += length(a)
    end

    return combined_indices
end

@doc raw"""
    minkowski_labels(A::AbstractVector{<:PointVector})

Calculates the labels of points in the Minkowski sum M of the point sets in A.
The points in M can be provided as optional argument if calculated beforehand.
"""
function minkowski_labels(
        A::AbstractVector{<:PointVector}...; 
        M=nothing,
        combined_indices=nothing
    )
    if M == nothing
        M = minkowski_sum(A...)
    end
    flatA = vcat(A...)

    #combined_indices = reduce(vcat, map(splat(fill), length.(A)|>enumerate))
    # combined_indices replaces every point in each A[i] with the index
    # this point would have if all A[i] were concatenated.
    if combined_indices == nothing
        combined_indices = points_to_total_indices(A...)
    end

    return Dict(map(Iterators.product(combined_indices...)) do v
        m = sum(flatA[collect(v)])

        v => findfirst(==(m), M)
    end)
end

@doc raw"""
    cayley_embedding(A::AbstractVector{<:PointVector}...)

For point sets $A_1, A_2, \dots, A_n\subset \mathbb{R}^d$, the Cayley embedding
is defined as $\mathcal{C}(A_1,\dots,A_n) = \bigcup_{j=1}^n \{ a\times e_j \mid a\in A_i \} \subset \mathbb{R}^d\times\mathbb{R}^n$.
"""
function cayley_embedding(A::AbstractVector{<:PointVector}...)
    n = length(A)

    return vcat(
        [vcat.(a, Ref([x==i for x in 1:n]))
         for (i,a) in enumerate(A)]...
    )
end

function inverse_cayley_embedding(C::AbstractVector{<:PointVector}, n::Int)
    A_with_j = map(x->(x[1:end-n], findfirst(==(1), x[end+1-n:end])), C)

    return map.(Ref(first), filter(x->last(x)==j, A_with_j) for j in 1:n)
end

@doc raw"""
    minkowski_projection(IncidenceMatrix, S::SubdivisionOfPoints, n::Int)

For a subdivision `S` of the Cayley embedding of `n` point sets, 
calculates the corresponding mixed subdivision of the Minkowski sum of the point sets
as incidence matrix.
"""
function minkowski_projection(::Type{IncidenceMatrix}, S::SubdivisionOfPoints, n::Int; M=nothing)
    return IncidenceMatrix.(
        minkowski_projection(points(S), maximal_cells(IncidenceMatrix, S), n; M=M)
    )
end

@doc raw"""
    minkowski_projection(IncidenceMatrix, S::SubdivisionOfPoints, n::Int)

For a subdivision `S` of the Cayley embedding of `n` point sets, 
calculates the corresponding mixed subdivision of the Minkowski sum of the point sets.
"""
function minkowski_projection(S::SubdivisionOfPoints, n::Int; M=nothing)
    if M == nothing
        C = points(S)
        A_with_j = map(C) do x 
            x[1:end-n], findfirst(==(1), x[end+1-n:end])
        end

        A = map(1:n) do j
            first.(filter(x->last(x)==j, A_with_j))
        end
        
        M = minkowski_sum(A...)
    end

    return subdivision_of_points(
             reduce(hcat, M) |> transpose, 
             minkowski_projection(IncidenceMatrix, S, n; M=M)
           )
end

@doc raw"""
    minkowski_projection(C::AbstractVector{<:PointVector}, cells::IncidenceMatrix, n::Int)

For a subdivision given by `cells` of the Cayley embedding `C` of `n` point sets,
calculates the corresponding mixed subdivision of the Minkowski sum.
"""
function minkowski_projection(
    C::AbstractVector{<:PointVector},
    cells::IncidenceMatrix,
    n::Int;
    M=nothing,
    labels=nothing
)
    A_with_j = map(C) do x 
        x[1:end-n], findfirst(==(1), x[end+1-n:end])
    end

    to_point_set = last.(A_with_j)
    A = map.(Ref(first), filter(x->last(x)==j, A_with_j) for j in 1:n)

    if M == nothing
        M = minkowski_sum(A...)
    end

    if labels == nothing
        labels = minkowski_labels(A...; M=M)
    end

    projected_incidences = map(1:nrows(cells)) do j
        separated_row = filter.([x->to_point_set[x]==i for i in 1:n], Ref(row(cells, j)))
        labels_of_cell = Iterators.product(separated_row...) |> collect |> vec

        get.(Ref(labels), labels_of_cell, 0)
    end

    return projected_incidences
end

@doc raw"""
    fine_mixed_subdivisions(as::Type{T} = SubdivisionOfPoints, G::Graph{Directed})
"""
function fine_mixed_subdivisions(as::Type{T}, G::Graph{Directed}; project_full=true) where {
  T<:Union{IncidenceMatrix, SubdivisionOfPoints, PolyhedralComplex}
}
  A = vertices_of_newton_polytope(G)

  if as == IncidenceMatrix
      return fine_mixed_subdivisions(IncidenceMatrix, A...)
  else 
      return fine_mixed_subdivisions(as, A...; project_full=project_full)
  end
end
fine_mixed_subdivisions(
    G::Graph{Directed}; 
    project_full=true
) = fine_mixed_subdivisions(SubdivisionOfPoints, G; project_full=project_full)

@doc raw"""
    fine_mixed_subdivisions(A::AbstractVector{<:PointVector}...)

Gives a list of fine mixed subdivisions of the Minkowski sum of `A[1]`, `A[2]`, ..., `A[n]`
together with the Minkowski sum itself.
"""
function fine_mixed_subdivisions(A::AbstractVector{<:PointVector}...)
    d = length(A |> first |> first)
    n = length(A)

    combined_indices = points_to_total_indices(A...)

    # Calculate which index belongs to which point set A_i
    M = minkowski_sum(A...)
    labels = minkowski_labels(A...; M=M, combined_indices=combined_indices)
    C = cayley_embedding(A...)

    CP = convex_hull(C)
    Tri = CP.pm_polytope |> Polymake.polytope.project_full |> polyhedron |> all_triangulations

    subdivisions = map(Tri) do T
      minkowski_projection(C, IncidenceMatrix(T), n; M=M, labels=labels)
    end

    return M, subdivisions
end

@doc raw"""
    fine_mixed_subdivisions(IncidenceMatrix, A::AbstractVector{<:PointVector}...)

Gives a list of incidence matrices for fine mixed subdivisions of the Minkowski sum 
of `A[1]`, `A[2]`, ..., `A[n]` together with the Minkowski sum itself.
"""
function fine_mixed_subdivisions(
    ::Type{IncidenceMatrix}, 
    A::AbstractVector{<:PointVector}...
)
    M, subdivisions = fine_mixed_subdivisions(A...)

    incidences = IncidenceMatrix.(subdivisions)
    for I in incidences
        resize!(I, nrows(I), length(M))
    end
    return M, incidences
end

@doc raw"""
    fine_mixed_subdivisions(SubdivisionOfPoints, A::AbstractVector{<:PointVector}...)

Gives a list of fine mixed subdivisions of the Minkowski sum of `A[1]`, `A[2]`, ..., `A[n]` 
as `SubdivisionOfPoints`.
"""
function fine_mixed_subdivisions(
    ::Type{SubdivisionOfPoints},
    A::AbstractVector{<:PointVector}...;
    project_full=true
)
    M, subdivisions = fine_mixed_subdivisions(IncidenceMatrix, A...)
    d = length(M[1])
    Mmat = reduce(hcat, M) |> transpose

    if project_full
      return subdivision_of_points.(Ref(project_matrix(Mmat)), subdivisions)
    else
      return subdivision_of_points.(Ref(Mmat), subdivisions)
    end
end

@doc raw"""
    fine_mixed_subdivisions(PolyhedralComplex, A::AbstractVector{<:PointVector}...)

Gives a list of fine mixed subdivisions of the Minkowski sum of `A[1]`, `A[2]`, ..., `A[n]` 
as `PolyhedralComplex`.
"""
function fine_mixed_subdivisions(
  ::Type{PolyhedralComplex},
  A::AbstractVector{<:PointVector}...;
  project_full=true
)
  M, subdivisions = fine_mixed_subdivisions(IncidenceMatrix, A...)

  Mmat = reduce(hcat, M) |> transpose
  d = length(M[1])

  if project_full
      return polyhedral_complex.(subdivisions, Ref(project_matrix(Mmat)))
  else 
      return polyhedral_complex.(subdivisions, Ref(Mmat))
  end
end

function polytrope_face_figures(::Type{IncidenceMatrix}, G::Graph{Directed})
  A = vertices_of_newton_polytope(G)

  M,T = fine_mixed_subdivisions(A...)

  polytrope_filters = filter.(Ref(cell->1 in cell), T) .|> IncidenceMatrix
  for I in polytrope_filters
      resize!(I, nrows(I), length(M))
  end

  return M, polytrope_filters
end
function polytrope_face_figures(G; project_full=true)
    return polytrope_face_figures(SubdivisionOfPoints, G; project_full=project_full)
end

function polytrope_face_figures(
  ::Type{SubdivisionOfPoints}, G::Graph{Directed}; project_full=true
)
  M, polytrope_filters = polytrope_face_figures(IncidenceMatrix, G)
  k = parent(M[1][1])
  d = length(M[1])

  Mmat = matrix(k, reduce(hcat, M) |> transpose)
  if project_full
      return subdivision_of_points.(Ref(project_matrix(Mmat)), polytrope_filters)
  else 
      return subdivision_of_points.(Ref(Mmat), polytrope_filters)
  end
end

function polytrope_face_figures(
  ::Type{PolyhedralComplex}, G::Graph{Directed}; project_full=true
)
  M, polytrope_filters = polytrope_face_figures(IncidenceMatrix, G)
  k = parent(M[1][1])
  d = length(M[1])

  Mmat = matrix(k, reduce(hcat, M) |> transpose)
  if project_full
      return polyhedral_complex.(polytrope_filters, Ref(project_matrix(Mmat)))
  else 
      return polyhedral_complex.(polytrope_filters, Ref(Mmat))
  end
end

function project_matrix(M::MatElem)
  d = size(M)[2]
  k = base_ring(M)

  p = sparse_matrix()
  for i in 2:d
      push!(p, sparse_row(k, [(i,k(1))]))
  end
  p = matrix(p) |> transpose
  return M*p
end

function secondary_fan_of_polytropes(G::Graph{Directed})
    A = cayley_embedding_of_dual_vertices(G)

    return secondary_fan(A)
end

function secondary_fan(V::AbstractVector{<:PointVector})
    k = parent(V[1][1])
    Vmat = matrix(k, reduce(hcat, V) |> transpose)

    return secondary_fan(Vmat)
end

function secondary_fan(V::MatElem)
    Polymake.Shell.tmp = convert(Polymake.PolymakeType, V)
    Polymake.shell_execute(raw"""$tmp = new polytope::PointConfiguration(POINTS=>$tmp);""")
    sF = Polymake.Shell.tmp |> Polymake.fan.secondary_fan |> polyhedral_fan
    Polymake.shell_execute(raw"""undef($tmp);""")

    return sF
end
