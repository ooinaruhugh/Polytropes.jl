using Oscar

function kleene_polynomials(T::Type{<:Vector}, G::Graph{Directed}) # for acyclic graphs
    I = indices(G)
    N = length(I)
    polynomials = Dict{Edge,Vector{Vector{Int}}}() # Save them as exponent vectors first, this is more economic for calculations

    for level in 1:nv(G)
        for j in 1:level-1
            # Translate edge into (exponent vector of) monomial
            c_ij = [zeros(Int, N)]
            c_ij[1][I[Edge(j,level)]] = 1 # this is just the variable `e_ij`
            
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
                isempty(c_kj) || append(c_kj,tmp)
            end

            setindex!(polynomials, c_ij, Edge(j,level))
        end
    end

    R = edge_ring(G)
    return get.(Ref(polynomials), edges_by_target(G), missing)
end

kleene_polynomials(T::Type{<:MPolyRing}, G::Graph{Directed}) = kleene_polynomials(Vector, G) .|> (x -> edge_ring(G)(ones(Int, length(x)), x))
kleene_polynomials(G::Graph{Directed}) = kleene_polynomials(MPolyRing, G)
end

function product_of_kleene_polynomials(G::Graph{Directed})
    return [filter(!ismonomial, kleene_polynomials(G)) |> prod]
end