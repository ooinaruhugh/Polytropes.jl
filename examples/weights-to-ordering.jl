using Oscar
using Polytropes

G = complete_dag(4)
R = edge_ring(G)
o = default_ordering(R)
I = toric_ideal(G)

GF = groebner_fan(I)

w = -Vector.(maximal_cones(GF) .|> rays_modulo_lineality .|> first .|> sum)
Q = weighted_digraph_polyhedron.(Ref(G), w)

Vector{Int}(v::Vector{QQFieldElem}) = Int.(v./gcd(v))
oW = weight_ordering.(Vector{Int}.(w),Ref(o))

F = kleene_polynomials(G)
inF = oW .|> (w->leading_term.(F; ordering=w)) .|> unique

