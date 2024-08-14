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

function action_on_edge_weights(G::Graph, s::PermGroupElem)
  N = n_edges(G)
  permuted_edges = map(edges(G)) do e
    Edge(e |> src |> s, e |> dst |> s)
  end |> enumerate |> collect |> x->sort(x; by=last)

  return first.(permuted_edges)
end

n = 4
function count_combinatorial_types(n)
  GG = transitively_closed_acyclic_graphs(n)
  II = all_pairs_shortest_path_ideal.(GG)

  gfan_with_data = map(zip(II,GG) |> collect |> filter(!iszero∘first)) do (I,G)
    I,G,groebner_fan(I)
  end

  count = length(GG) - length(gfan_with_data)

  count += map(gfan_with_data) do (_,G,F)
    int_rays = Polytropes.interior_points_of_cones(F)
    Aut = automorphism_group(G)

    if is_trivial(Aut)
      length(int_rays)
    else
      map(S) do s
        map(int_rays) do v
          map(maximal_cones(F)) do C
            v[action_on_edge_weights(G,s)] in C
          end
        end
      end |> splat(zip) .|> sum |> unique |> length
    end
  end |> sum

  return count
end

#C_subdivisions = map(gfan_with_data) do (I,G,gfan)
#  W = Polytropes.interior_points_of_cones(gfan)
#  Q = weighted_digraph_polyhedron.(Ref(G), W)
#
#  A = Polytropes.vertices_of_newton_polytope(G)
#  C = Polytropes.cayley_embedding(A...)
#  H = extend_weight_to_height.(Ref(G), W)
#  subdivision_of_points.(Ref(C), H)
#end
#
