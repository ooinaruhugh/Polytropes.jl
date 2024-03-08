module Transitivity

using Graphs
const Graphs_jl = Graphs
using Oscar

function to_graphs_graph(G::Oscar.Graph{Directed})
    out = DiGraph(nvertices(G))
    for e in Oscar.edges(G)
        Graphs_jl.add_edge!(out, Oscar.src(e), Oscar.dst(e))
    end
    
    return out
end

function to_oscar_graph(G::DiGraph)
    out = Oscar.Graph{Directed}(Graphs_jl.nv(G))
    for e in Graphs_jl.edges(G)
        Oscar.add_edge!(out, Graphs_jl.src(e), Graphs_jl.dst(e))
    end
    
    return out
end

function transitive_closure(G::Oscar.Graph{Directed})
    H = to_graphs_graph(G)
    transitiveclosure!(H)

    return to_oscar_graph(H)
end

function transitive_reduction(G::Oscar.Graph{Directed})
    H = to_graphs_graph(G)
    out = transitivereduction(H)

    return to_oscar_graph(out)
end

export transitive_closure
export transitive_reduction

end