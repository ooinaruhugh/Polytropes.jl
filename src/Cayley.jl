using Oscar

function cayley_embedding_of_dual_vertices(G::Graph{Directed})
    C = Vector{Vector{Int}}[]
    n = n_vertices(G)

    for v in 1:n
        Cv = Vector{Int}[]
        for u in [v,outneighbors(G, v)...]
            c = zeros(Int, n)
            c[u] = 1

            push!(Cv,c)
        end

        push!(C,Cv)
    end

    return cayley_embedding(C...)
end

function cayley_embedding(A::AbstractVector{<:AbstractVector{T}}...) where {T}
    C = Vector{Int}[]
    d = length(A |> first |> first)
    n = length(A)
    for (i,a) in enumerate(A)
        for p in a
            c = [p...,zeros(T, n)...]
            c[d+i] = 1

            push!(C,c)
        end
    end

    return C
end

function minkowski_sum(A::AbstractVector{<:AbstractVector{T}}...) where {T}
    return vec([sum(v) for v in Iterators.product(A...)])
end

function mixed_subdivisions(G::Graph{Directed})
  n = n_vertices(G)
  f = x->[(x,1)]
  
  A = [
       Vector.(sparse_row.(Ref(ZZ), f.([v,outneighbors(G,v)...])) ,n)
    for v in 1:n
  ]

  return mixed_subdivisions(A...)
end

function mixed_subdivisions(A::AbstractVector{<:AbstractVector{T}}...) where {T}
    d = length(A |> first |> first)
    n = length(A)

    _indices = [
        Ref(l).+(1:length(a)) 
        for (a,l) in Iterators.zip(A,
            [0,(sum(length.(A[1:i-1])) for i in 2:n)...]
        )
    ]
    minkowski_lookup = Dict(
        t => i for (i,t) in enumerate(Iterators.product(_indices...))
    )
            
    M = minkowski_sum(A...)
    C = cayley_embedding(A...)
    cayley_lookup = Dict(i => findfirst(x->x==1, c[d+1:d+n]) for (i,c) in enumerate(C))

    CP = convex_hull(C)
    Tri::Vector{Vector{Vector{Int}}} = CP.pm_polytope |> Polymake.polytope.project_full |> polyhedron |> all_triangulations

    subdivisions = Vector{Vector{Vector{Int}}}[]
    for t in Tri
        subdivision = [
            [
                minkowski_lookup[I]
                for I in vec(collect(Iterators.product((filter(x->cayley_lookup[x] == i, Δ) for i in 1:n)...)))
            ] for Δ in t
        ]

        push!(subdivisions,[subdivision])
    end

    return M, subdivisions
end

function polytrope_duals(G::Graph{Directed})
  M,T = mixed_subdivisions(G)

  return mixed_subdivisions_as_complexes(M, unique(filter.(x->1∈x, t) for t in T))
end

mixed_subdivisions_as_complexes(M, subdivisions) = [
  polyhedral_complex(IncidenceMatrix(subdivision...), M) 
  for subdivision in subdivisions
]

mixed_subdivisions_as_complexes(
  A::AbstractVector{<:AbstractVector{T}}...
) where {T} = mixed_subdivisions_as_complexes(mixed_subdivisions(A...))
# return [
#     polyhedral_complex(IncidenceMatrix(subdivision...), M).pm_complex |> Polymake.polytope.project_full |> polyhedral_complex
#     for subdivision in subdivisions
# ]
