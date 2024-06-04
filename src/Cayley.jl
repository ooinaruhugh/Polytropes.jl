using Oscar

function vertices_of_newton_polytope(G::Graph{Directed})
    C = Vector{PointVector}[]
    n = n_vertices(G)

    for v in 1:n
        Cv = [
            point_vector([i==u for i in  1:n])
            for u in [v,outneighbors(G, v)...]
        ]

        push!(C,Cv)
    end

    return C
end

function cayley_embedding_of_dual_vertices(G::Graph{Directed})
    C = vertices_of_newton_polytope(G)
    return cayley_embedding(C...)
end

function cayley_embedding(A::AbstractVector{<:PointVector}...)
    n = length(A)

    return vcat(
        [vcat.(a, Ref([x==i for x in 1:n]))
         for (i,a) in enumerate(A)]...
    )
end

function coherent_minkowski_indices(
        A::AbstractVector{<:AbstractVector{<:PointVector}}, 
        M::AbstractVector{<:PointVector}
    )
    n = length(A)
    d = length(M[1])

    _indices = [
        Ref(l).+(1:length(a)) 
        for (a,l) in Iterators.zip(A,
            [0,(sum(length.(A[1:i-1])) for i in 2:n)...]
        )
    ]

    C = vcat(A...)
    lookup = Dict(
                  t => findfirst(==(sum(C[collect(t)])), M)
                  for t in Iterators.product(_indices...)
                 )

    return lookup
end

function cayley_to_minkowski_subdivision(
        C::AbstractVector{<:PointVector},
        cells::IncidenceMatrix
    )
    A = minkowski_projection(C)
    M = minkowski_sum(A...)

    labels = coherent_minkowski_indices(A, M)

    _indices = Dict(vcat(
        [Ref(l).+(1:length(a)) .=> Ref(i)
        for ((i,a),l) in Iterators.zip(enumerate(A),
            [0,(sum(length.(A[1:i-1])) for i in 2:n)...]
        )
       ]...)...)

    return _cayley_to_minkowski_subdivision(M, cells, labels, _indices)
end

function _cayley_to_minkowski_subdivision(
        M::AbstractVector{<:PointVector},
        cells::AbstractVector{<:Vector},
        cayley_to_minkowski::Dict,
        which_point_set::Dict
    )
    incidence = IncidenceMatrix([
             vec([lookup[v] for v in Iterators.product(t...)]) for t in separated_T
            ])
end

function minkowski_sum(A::AbstractVector{<:PointVector}...)
    return vec([sum(v) for v in Iterators.product(A...)]) |> unique
end

function minkowski_projection(C::AbstractVector{<:PointVector}, n::Int)
  d = length(C[1]) - n
  A = [ 
       [ C[j][1:d] for j in 1:length(C) if findfirst(==(1), C[j][end-n+1:end]) == i] 
        for i in 1:n 
      ]

  return A
end

function mixed_subdivisions(G::Graph{Directed})
  A = vertices_of_newton_polytope(G)
  return mixed_subdivisions(A...)
end

function mixed_subdivisions(::Type{IncidenceMatrix}, G::Graph{Directed})
  A = vertices_of_newton_polytope(G)
  return mixed_subdivisions(IncidenceMatrix, A...)
end

function mixed_subdivisions(A::AbstractVector{<:PointVector}...)
    d = length(A |> first |> first)
    n = length(A)

    _indices = Dict(vcat(
        [Ref(l).+(1:length(a)) .=> Ref(i)
        for ((i,a),l) in Iterators.zip(enumerate(A),
            [0,(sum(length.(A[1:i-1])) for i in 2:n)...]
        )
       ]...)...)

    # Calculate which index belongs to which point set A_i
    M = minkowski_sum(A...)
    minkowski_lookup = coherent_minkowski_indices(collect(A),M)
    C = cayley_embedding(A...)

    CP = convex_hull(C)
    Tri::Vector{Vector{Vector{Int}}} = CP.pm_polytope |> Polymake.polytope.project_full |> polyhedron |> all_triangulations

    subdivisions = Vector{<:Vector}[]
    for T in Tri
      separated_T = [filter.([v->_indices[v]==i for i in 1:n], Ref(t)) for t in T]

      push!(subdivisions, 
            [
             vec([minkowski_lookup[v] for v in Iterators.product(t...)]) for t in separated_T
            ])
    end

    return M, subdivisions
end

function mixed_subdivisions(::Type{IncidenceMatrix}, A::AbstractVector{<:PointVector}...)
    M, subdivisions = mixed_subdivisions(A...)

    incidences = IncidenceMatrix.(subdivisions)
    for I in incidences
        resize!(I, nrows(I), length(M))
    end
    return M, IncidenceMatrix.(subdivisions)
end

function polytrope_duals(G::Graph{Directed})
  M,T = mixed_subdivisions(G)
  incidences = IncidenceMatrix.([filter(x->1∈x, t) for t in T])
  for I in incidences
      resize!(I, nrows(I), length(M))
  end

  return mixed_subdivisions_as_complexes(M, incidences)
end

function mixed_subdivisions_as_complexes(M::AbstractVector{<:PointVector}, subdivisions; project_full=true) 
  Mmat = reduce(hcat, M) |> transpose

  complexes = polyhedral_complex.(subdivisions, Mmat)
  if project_full
      output = PolyhedralComplex[]
      for P in complexes
          maximal_polyhedra(P)
          push!(output, P.pm_complex |> Polymake.polytope.project_full |> polyhedral_complex)
      end

      return output
  else return complexes
  #return subdivision_of_points.(Ref(Mmat), subdivisions)
end

mixed_subdivisions_as_complexes(
  A::AbstractVector{<:AbstractVector{T}}...;
  project_full=true
) where {T} = mixed_subdivisions_as_complexes(mixed_subdivisions(IncidenceMatrix, A...); project_full=true)
# return [
#     polyhedral_complex(IncidenceMatrix(subdivision...), M).pm_complex |> Polymake.polytope.project_full |> polyhedral_complex
#     for subdivision in subdivisions
# ]
