using Oscar
using Polytropes

# Example 1: Enumerate the combinatorial types for a complete dag on 4 nodes
G = complete_dag(4)
K, F = kleene_graph(G)

# This should give 9 solutions
S = enumerate_satisfying_assignments(K,F)

# Example 2: Given a solution for a supergraph, we can recover all solutions without running the above computation again
# For this
#   - we need to remove the edges we don't use
#   - prune solutions that use those edges
#   - take care of graph automorphisms that our subgraph might admit

# Here, we do this for the diamond dag sitting inside the complete dag.
# For this, we need to remove the edge 2->3

v = 4
Kop = opposite_graph(K)
KopT = transitive_closure(Kop)

S1 = S .|> (s->filter(x->x!=v, s)) |> (s->filter(x->intersect(x, outneighbors(KopT, v)) |> isempty, s))

M = monomials_of_kleene_polynomials(G)
S1 .|> s->M[s]