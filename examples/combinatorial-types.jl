using Oscar
using Polytropes

S = complete_dag(4) |> enumerate_satisfying_assignments

# Example 1: Enumerate the combinatorial types for a complete dag on 4 nodes
G = complete_dag(4)
K, F = kleene_graph(G)

# This should give 9 solutions
S = enumerate_satisfying_assignments(K,F) .|> sort!

# Example 2: Given a solution for a supergraph, we can recover all solutions without running the above computation again
# For this
#   - we need to remove the edges we don't use
#   - prune solutions that use those edges
#   - take care of graph automorphisms that our subgraph might admit

# Here, we do this for the diamond dag sitting inside the complete dag.
# For this, we need to remove the edge 2->3 (which has number 4 in this specific Kleene graph)

v = 4
Kop = opposite_graph(K)
KopT = transitive_closure(Kop)

contains_v = [v, outneighbors(KopT, v)...]
intersects = s -> (x -> intersect(x,s) |> !isempty)

prune_v = s -> filter(x->x!=v, s)

S1 = S .|> prune_v |> (s->filter(!intersects(contains_v), s))

M = monomials_of_kleene_polynomials(G)
S1 .|> s->M[s]
# This gives the three solutions you would expect!

R = edge_ring(G)
M = monomials_of_kleene_polynomials(G)
a = findfirst(x->x == M[v], gens(R)) |> i -> edges_by_target(G)[i]

H = graph_from_edges(Directed, [e for e in edges(G) if e != a])
P = edge_ring(H)
index_in_subgraph = e -> findfirst(x->x==e, edges_by_target(H))
gen_in_subring = i -> gens(P)[i]

im_π = [has_edge(H, e) ? e |> index_in_subgraph |> gen_in_subring : 0 for e in edges_by_target(G)]
Π = hom(R, P, im_π)

subindices = [i for (i,m) in enumerate(Π.(M)) if m!=0]

σ = automorphism_group_generators(H) |> first
            

im = [Edge(e |> src |> σ, e |> dst |> σ) for e in edges(H)] .|> 
         e -> findfirst(y->y==e, edges_by_target(H)) .|>
         i -> gens(P)[i]

ϕ = hom(P,P, im)

filter_0 = v -> filter(x->x!=0, v)
α = Dict(zip(subindices, ϕ.(Π.(M)) |> filter_0 .|> m -> findfirst(x->x==m, Π.(M))))

S2 = S1 .|> s->[α[i] for i in s]

# TODO: Graph automorphisms and action on K