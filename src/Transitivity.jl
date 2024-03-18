module Transitivity

using Graphs
const _Graphs = Graphs
using Oscar

function to_graphs_graph(G::Oscar.Graph{Directed})
    out = DiGraph(n_vertices(G))
    for e in Oscar.edges(G)
        _Graphs.add_edge!(out, Oscar.src(e), Oscar.dst(e))
    end
    
    return out
end

function to_oscar_graph(G::DiGraph)
    out = Oscar.Graph{Directed}(_Graphs.nv(G))
    for e in _Graphs.edges(G)
        Oscar.add_edge!(out, _Graphs.src(e), _Graphs.dst(e))
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