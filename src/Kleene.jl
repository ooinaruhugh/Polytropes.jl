using Oscar

function kleene_polynomials(G::Graph{Directed})
    I = indices(G)
    polynomials = Dict{Edge,Vector{Vector{Int}}}()

    for level in 1:nv(G)
        for j in 1:level 
            if j != level
                c_ij = [[(Edge(j,level) == I(e)) ? 1 : 0 for e in 1:length(I)]]
                
                for k in inneighbors(G, level)
                    c_ik = get(polynomials, Edge(j,k), [])
                    e_k = get(I, Edge(k,level), 0)
                    
                    tmp = [[i!=e_k ? a : a+1 for (i,a) in enumerate(m)] for m in c_ik]
                    isempty(c_ik) || append!(c_ij,tmp)
                end
                # for k in outneighbors(G, j)
                #     c_kj = get(polynomials, Edge(k,level), [])
                #     e_k = get(I, Edge(j,level), 0)
                    
                #     tmp = [[i!=e_k ? a : 1 for (i,a) in enumerate(m)] for m in c_kj]
                #     isempty(c_kj) || append(c_kj,tmp)
                # end

                setindex!(polynomials, c_ij, Edge(j,level))
            end 
        end
    end

    ctx = MPolyBuildCtx(edge_ring(G))
    return map(x -> (x .|> (term -> push_term!(ctx, 1, term)), finish(ctx)) |> last,
        polynomials |> values)
        
end

function product_of_kleene_polynomials(G::Graph{Directed})
    return [filter(!ismonomial, kleene_polynomials(G)) |> prod]
end