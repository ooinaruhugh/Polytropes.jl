using Oscar
# import Graphs: DiGraph, transitiveclosure, add_edge!, nv
# using MetaGraphs

function unit_vector(T::Type, N, i)
    v = zeros(T, N)
    v[i] = 1
    return v
end

function kleene_polynomials(::Type{<:MPolyRing}, G::Graph{Directed}) # for acyclic graphs
    I = indices(G)
    N = length(I)
    polynomials = Dict{Edge,Vector{Vector{Int}}}() # Save them as exponent vectors first, this is more economic for calculations
    transitive_edges = []

    for level in 1:nv(G)
        for j in 1:level
            c_ij = []
            has_edge(G, j, level) && push!(c_ij, unit_vector(Int, N, I[Edge(j,level)]))
            
            for k in inneighbors(G, level) # costruct summands of ij-th polynomial as product of ik-th edge times kj-th polynomial
                c_ik = get(polynomials, Edge(j,k), []) 
                e_k = get(I, Edge(k,level), 0)
                
                tmp = [[i!=e_k ? a : a+1 for (i,a) in enumerate(m)] for m in c_ik]
                isempty(c_ik) || append!(c_ij,tmp)
            end
            # The following is apparently not necessary
            for k in outneighbors(G, j)
                c_kj = get(polynomials, Edge(k,level), [])
                e_k = get(I, Edge(j,level), 0)
                
                tmp = [[i!=e_k ? a : 1 for (i,a) in enumerate(m)] for m in c_kj]
                isempty(c_kj) || append!(c_kj,tmp)
            end

            setindex!(polynomials, c_ij, Edge(j,level))
            !isempty(c_ij) && push!(transitive_edges, Edge(j,level))
        end
    end

    R = edge_ring(G)
    return get.(Ref(polynomials), transitive_edges, missing) .|> (x -> R(ones(Int, length(x)), x))
end

# We need to convert from exponent vectors to polynomials and then back so the ordering stays stable
kleene_polynomials(::Type{<:Vector}, G::Graph{Directed}) = kleene_polynomials(G) .|> exponents .|> collect
kleene_polynomials(G::Graph{Directed}) = kleene_polynomials(MPolyRing, G)

function kleene_graph(G::Graph{Directed})
    N = nedges(G)

    F = kleene_polynomials(Vector, G)
    V = reduce(vcat, F)
    groups = sort(F;lt=!isless,by=length) .|> (f -> [findfirst(x->x==m,V) for m in f])

    # Find the index of the first and last edge in every path. The variables are sorted nicely in acyclic graphs, so they're just the first and last non-zero exponent.
    k = V .|> y->[findfirst(x->x==1, y), findlast(x->x==1, y)]
    factored_paths::Vector{Vector{Vector{Int}}} = [
        [x - unit_vector(Int, N, k[i] |> first), x - unit_vector(Int, N, k[i] |> last)]
        for (i,x) in enumerate(V)
    ]

    # `factored_paths` is a Vector of pairs of monomials/exponent vectors, so we apply `findfirst` to both components z of every pair x
    EK = factored_paths .|> x->(x.|> z -> findfirst(y->y==z,V))

    K = Graph{Directed}(length(V))
    for (i,x) in enumerate(EK)
        first(x) |> isnothing || add_edge!(K, i, first(x))
        last(x) |> isnothing || add_edge!(K,i,last(x))
    end
    return K, groups
end

function build_reverse_lookup(F::Vector{Vector{Int}})
    lookup = Dict{Int, Vector{Int}}()
    for f in F
        for m in f
            lookup[m] = f
        end
    end

    return lookup
end

function enumerate_satisfying_assignments(K::Graph{Directed}, F::Vector{Vector{Int}}, lookup::Dict)
    f = first(F)
    solutions = []

    for m in f
        solution = [m, neighbors(K,m)...]

        assigned = [get(lookup, x, nothing) for x in solution]
        unassigned = filter(x->x ∉ assigned, F)

        if unassigned |> !is_empty
            subsolutions = enumerate_satisfying_assignments(K, unassigned, lookup)

            for subsolution in subsolutions
                push!(solutions, union(solution, subsolution))
            end
        else
            push!(solutions, solution)
        end
    end
    
    return solutions
end

function enumerate_satisfying_assignments(K::Graph{Directed}, F::Vector{Vector{Int}})
    lookup = build_reverse_lookup(F)
    KT = transitive_closure(K)

    return enumerate_satisfying_assignments(KT, F, lookup)    
end

function product_of_kleene_polynomials(G::Graph{Directed})
    return filter(!ismonomial, kleene_polynomials(G)) |> prod
end