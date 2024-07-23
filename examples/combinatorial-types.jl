using Oscar
using Polytropes

function extend_weight_to_height(G::Graph, w::T) where {U, T<:AbstractVector{U}}
  h = U[]

  ew = zip(edges(G), w) |> collect

  for v in 1:n_vertices(G)
      push!(h, 0)
      for (e,w) in ew
        if src(e) == v
            push!(h,w)
        end
      end
  end

  return h
end

n = 4
GG = transitively_closed_acyclic_graphs(n)
II = all_pairs_shortest_path_ideal.(GG)

gfan_with_data = map(zip(II,GG) |> collect |> filter(!iszero∘first)) do (I,G)
  I,G,groebner_fan(I)
end

C_subdivisions = map(gfan_with_data) do (I,G,gfan)
  W = Polytropes.interior_points_of_cones(gfan)
  Q = weighted_digraph_polyhedron.(Ref(G), W)

  A = Polytropes.vertices_of_newton_polytope(G)
  C = Polytropes.cayley_embedding(A...)
  H = extend_weight_to_height.(Ref(G), W)
  subdivision_of_points.(Ref(C), H)
end

