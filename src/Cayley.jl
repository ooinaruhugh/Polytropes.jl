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

function minkowski_sum(A::AbstractVector{<:PointVector}...)
    return vec([sum(v) for v in Iterators.product(A...)]) |> unique
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
    return M, IncidenceMatrix.(subdivisions)
end

function polytrope_duals(G::Graph{Directed})
  M,T = mixed_subdivisions(G)

  return mixed_subdivisions_as_complexes(M, IncidenceMatrix.([filter(x->1∈x, t) for t in T]))
end

function mixed_subdivisions_as_complexes(M, subdivisions) 
  Mmat = reduce(hcat, M) |> transpose
  return subdivision_of_points.(Ref(Mmat), subdivisions)
end

mixed_subdivisions_as_complexes(
  A::AbstractVector{<:AbstractVector{T}}...
) where {T} = mixed_subdivisions_as_complexes(mixed_subdivisions(A...))
# return [
#     polyhedral_complex(IncidenceMatrix(subdivision...), M).pm_complex |> Polymake.polytope.project_full |> polyhedral_complex
#     for subdivision in subdivisions
# ]
