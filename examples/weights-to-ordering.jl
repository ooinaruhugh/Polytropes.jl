using Oscar
using Polytropes

G = complete_dag(3)
R = edge_ring(G)
I = toric_ideal(G)

# w = [-1,-2,-3]
w = [-1,-3,-2]

o = lex(R)
oW = weight_ordering(w,o)

leading_ideal(I;ordering=oW)