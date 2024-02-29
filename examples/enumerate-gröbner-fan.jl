using Oscar
using Polytropes

Vector{ZZRingElem}(v::RayVector{QQFieldElem}) = ZZ.(convert(Vector, v))
Vector{Int}(v::RayVector{QQFieldElem}) = Int.(convert(Vector, v))

n = 4

G = complete_dag(n)
I = toric_ideal(G)

GF = groebner_fan(I)
w = GF |> maximal_cones .|> rays_modulo_lineality .|> first .|> sum
int_w = w .|> (x->(denominator.(x) |> lcm)*x) .|> Vector{Int}

Q = w .|> x -> weighted_digraph_polyhedron(G,x)

f = f_vector.(Q)

ins = int_w .|> (x->weight_ordering(x,lex(edge_ring(G)))) .|> x->leading_ideal(I;ordering=x)