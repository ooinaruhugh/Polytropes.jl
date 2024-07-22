function transitively_closed_acyclic_graphs(n::Int)
    output = Graph{Directed}[]
    seen = Set{Int}()

    queue = Graph{Directed}[complete_dag(n)]
    while !is_empty(queue)
        G = popfirst!(queue)
        push!(output, G)
        E = transitive_reduction(G) |> edges |> collect

        for e in E
            rem_edge!(G, src(e), dst(e))

            graph_hash = Polymake.graph.canonical_hash(Oscar.pm_object(G))
            if graph_hash ∉ seen
              if edges(G) |> !is_empty
                H = graph_from_edges(Directed, edges(G), n)
                push!(queue, H)
              end
            end
            push!(seen, graph_hash)

            add_edge!(G, src(e), dst(e))
        end
    end

    push!(output, Graph{Directed}(n))

    return output
end

function fan_modulo_lineality(P::PolyhedralFan)

end

function interior_points_of_cones(P::PolyhedralFan)
  if lineality_dim(P) > 0
      return maximal_cones(P) .|> rays_modulo_lineality .|> first .|> sum
  else
      return maximal_cones(P) .|> rays .|> sum
  end
end
