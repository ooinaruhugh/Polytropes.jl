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

    labels = Pair[]
    for v in Iterators.product(combined_indices...)
        m = sum(flatA[collect(v)])

        push!(labels, v => findfirst(==(m), M))
    end

    return Dict(labels)
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

    return map.(Ref(first), 
             filter(x->last(x)==j, A_with_j) for j in 1:n
            )
end

@doc raw"""
    minkowski_projection(IncidenceMatrix, S::SubdivisionOfPoints, n::Int)

For a subdivision `S` of the Cayley embedding of `n` point sets, 
calculates the corresponding mixed subdivision of the Minkowski sum of the point sets
as incidence matrix.
"""
function minkowski_projection(::Type{IncidenceMatrix}, S::SubdivisionOfPoints, n::Int; M=nothing)
    return minkowski_projection(points(S), maximal_cells(IncidenceMatrix, S), n; M=M)
end

@doc raw"""
    minkowski_projection(IncidenceMatrix, S::SubdivisionOfPoints, n::Int)

For a subdivision `S` of the Cayley embedding of `n` point sets, 
calculates the corresponding mixed subdivision of the Minkowski sum of the point sets.
"""
function minkowski_projection(S::SubdivisionOfPoints, n::Int; M=nothing)
    if M == nothing
        C = points(S)
        A_with_j = map(x->(x[1:end-n], findfirst(==(1), x[end+1-n:end])), C)

        A = map.(Ref(first), 
                 filter(x->last(x)==j, A_with_j) for j in 1:n
                )
        
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

    A_with_j = map(x->(x[1:end-n], findfirst(==(1), x[end+1-n:end])), C)

    to_point_set = last.(A_with_j)
    A = map.(Ref(first), 
             filter(x->last(x)==j, A_with_j) for j in 1:n
            )

    if M == nothing
        M = minkowski_sum(A...)
    end

    if labels == nothing
        labels = minkowski_labels(A...; M=M)
    end

    projected_incidences = Vector{Int}[]
    for j in 1:nrows(cells)
        separated_row = filter.([x->to_point_set[x]==i for i in 1:n], Ref(row(cells, j)))
        labels_of_cell = Iterators.product(separated_row...) |> collect |> vec
        push!(projected_incidences, get.(Ref(labels), labels_of_cell, 0))
    end

    return IncidenceMatrix(projected_incidences)
end

function fine_mixed_subdivisions(G::Graph{Directed})
  A = vertices_of_newton_polytope(G)
  return fine_mixed_subdivisions(A...)
end

@doc raw"""
    fine_mixed_subdivisions(IncidenceMatrix, G::Graph{Directed})
"""
function fine_mixed_subdivisions(::Type{IncidenceMatrix}, G::Graph{Directed})
  A = vertices_of_newton_polytope(G)
  return fine_mixed_subdivisions(IncidenceMatrix, A...)
end

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

    subdivisions = [
      minkowski_projection(C, IncidenceMatrix(T), n; M=M, labels=labels)
      for T in Tri
    ]

    return M, subdivisions
end

@doc raw"""
    fine_mixed_subdivisions(IncidenceMatrix, A::AbstractVector{<:PointVector}...)

Gives a list of `IncidenceMatrix` for fine mixed subdivisions of the Minkowski sum 
of `A[1]`, `A[2]`, ..., `A[n]` together with the Minkowski sum itself.
"""
function fine_mixed_subdivisions(::Type{IncidenceMatrix}, A::AbstractVector{<:PointVector}...)
    M, subdivisions = mixed_subdivisions(A...)

    incidences = IncidenceMatrix.(subdivisions)
    for I in incidences
        resize!(I, nrows(I), length(M))
    end
    return M, IncidenceMatrix.(subdivisions)
end

function polytrope_duals(G::Graph{Directed}; project_full=false)
  M,T = mixed_subdivisions(G)
  d = length(M[1])
  incidences = IncidenceMatrix.([filter(x->1∈x, t) for t in T])
  for I in incidences
      resize!(I, nrows(I), length(M))
  end

  Mmat = matrix(QQ, reduce(hcat, M) |> transpose)
  if project_full
      π = sparse_matrix(QQ)
      for i in 2:d
          push!(π, sparse_row(QQ, [(i,QQ(1))]))
      end
      π = matrix(π) |> transpose

      return polyhedral_complex.(incidences, Ref(Mmat*π))
  else 
      return polyhedral_complex.(incidences, Ref(Mmat))
  end
end

function mixed_subdivisions_as_complexes(M::AbstractVector{<:PointVector}, subdivisions; project_full=true) 
  Mmat = reduce(hcat, M) |> transpose

  complexes = polyhedral_complex.(subdivisions, Ref(Mmat))
  if project_full
      output = PolyhedralComplex[]
      for P in complexes
          maximal_polyhedra(P)
          push!(output, P.pm_complex |> Polymake.polytope.project_full |> polyhedral_complex)
      end

      return output
  else return complexes
  end
  #return subdivision_of_points.(Ref(Mmat), subdivisions)
end

mixed_subdivisions_as_complexes(
  A::AbstractVector{<:AbstractVector{T}}...;
  project_full=true
) where {T} = mixed_subdivisions_as_complexes(mixed_subdivisions(IncidenceMatrix, A...); project_full=true)
