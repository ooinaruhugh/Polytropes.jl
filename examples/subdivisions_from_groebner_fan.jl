using Oscar
using Polytropes

n = 4
G = complete_dag(4)

C = cayley_embedding_of_dual_vertices(G)
Cmat = matrix(QQ, reduce(hcat, C) |> transpose)

sF = secondary_fan_of_polytropes(G)
h = interior_points_of_cones(sF)

Q_from_secondary_fan = weighted_digraph_polyhedron.(Ref(G), height_to_weight.(Ref(G), h))
duals_from_secondary_fan = subdivision_of_points.(Ref(Cmat), h)

I = toric_ideal(G)
GF = groebner_fan(I)
w = interior_points_of_cones(GF)

Q_from_gfan = weighted_digraph_polyhedron.(Ref(G), w)
duals_from_gfan = subdivision_of_points.(Ref(Cmat), weight_to_height.(Ref(G), w))
