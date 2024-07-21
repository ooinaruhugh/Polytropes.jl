using Oscar
using Polytropes

G = complete_dag(4)
R = edge_ring(G)
o = default_ordering(R)
I = toric_ideal(G)
F = kleene_polynomials(G)

make_integral_coords(v::AbstractVector{QQFieldElem}) = Int.(v./gcd(v))

# Method 1: Using the Groebner fan to produce weights
## Creating the weighted digraph polyhedra
GF = groebner_fan(I)
GF_w = Vector.(maximal_cones(GF) .|> rays_modulo_lineality .|> first .|> sum)
GF_Q = weighted_digraph_polyhedron.(Ref(G), GF_w)

## Comparing with the satisfying assignments
oW = weight_ordering.(make_integral_coords.(GF_w),Ref(o))
inF = oW .|> (w->leading_term.(F; ordering=w)) 

# Method 2: Using the secondary fan of the Cayley embedding
C = Polytropes.cayley_embedding_of_dual_vertices(G)

Polymake.Shell.V = convert(Polymake.PolymakeType, QQ.(reduce(hcat,C)'))
Polymake.shell_execute(raw"""application "fan";""")
Polymake.shell_execute(raw"""$p = new PointConfiguration(POINTS=>$V);""")
sF = Polymake.fan.secondary_fan(Polymake.Shell.p) |> polyhedral_fan

indices = [2,3,4,6,7,9]
interior_pts = maximal_cones(sF) .|> rays_modulo_lineality .|> first .|> sum .|> Vector
interior_matrix = reduce(hcat, interior_pts) |> transpose

cayley_Q = 1:20 .|> (i -> weighted_digraph_polyhedron(G, interior_matrix[i, indices]))


