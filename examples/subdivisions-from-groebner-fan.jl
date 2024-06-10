using Oscar
using Polytropes

n = 4
G = complete_dag(n)

C = cayley_embedding_of_dual_vertices(G)
Cmat = matrix(QQ, reduce(hcat, C) |> transpose)

sF = secondary_fan(C)
h = Polytropes.interior_points_of_cones(sF)

Q_from_secondary_fan = weighted_digraph_polyhedron.(Ref(G), Polytropes.height_to_weight.(Ref(G), h))
duals_from_secondary_fan = map(subdivision_of_points.(Ref(Cmat), h)) do s
    Ms = minkowski_projection(s,n)
    Polytropes.embed_subdivision_in_full_simplex(Ms, n)
end

I = toric_ideal(G)
GF = groebner_fan(I)
w = Polytropes.interior_points_of_cones(GF)

Q_from_gfan = weighted_digraph_polyhedron.(Ref(G), w)
duals_from_gfan = map(subdivision_of_points.(Ref(Cmat), Polytropes.weight_to_height.(Ref(G), w))) do s
    Ms = minkowski_projection(s,n)
    Polytropes.embed_subdivision_in_full_simplex(Ms, n)
end

