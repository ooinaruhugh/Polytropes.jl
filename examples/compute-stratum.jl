using Oscar, Polytropes
using CodecXz

function graph_from_fundamental_polytope(A::T) where {T <: Union{Matrix,MatElem}}
  E = map(eachrow(A[2:end,2:end])) do v
    s = findfirst(==(-1), v)
    t = findfirst(==(1) , v) + 1
    s = isnothing(s) ? 1 : s+1

    Edge(s,t)
  end

  return graph_from_edges(Directed, E)
end

function compute_stratum(input_file::String)
  if endswith(input_file, ".xz")
    triangs = open(`xzcat $input_file`, "r") do io
      load(io)
    end
  else
    triangs = load(input_file)
  end
  output_file = rsplit(input_file, ".mrdi"; limit=2)[1] * ".fan.mrdi.xz"

  D = triangs.points |> n_rows

  G = graph_from_fundamental_polytope(triangs.points)

  Σ = map(triangs.triangulations) do t
    subdivision_of_points(triangs.points, t)
  end
  scones = secondary_cone.(Σ)

  if lineality_dim(scones[1]) < ambient_dim(scones[1])
    L, = lineality_space.(scones) |> unique
    L = reduce(hcat, L) |> transpose |> matrix
    _,rrefL = rref(L)
    i = map(eachrow(rrefL)) do v
      findfirst(!iszero,v)
    end

    red_rays = map(scones) do C
      R,_ = rays_modulo_lineality(C)
      rr = map(R) do r
        λ = r[i] |> collect
        q = λ * rrefL
        (r - q)[setdiff(1:D,i)]
      end

      reduce(hcat, rr) |> transpose |> matrix
    end

    red_scones = map(zip(scones,red_rays)) do (scone,r)
      polyhedral_fan(cones(scone|>polyhedral_fan),r; non_redundant=true) |> maximal_cones |> first
    end

    red_sfan = polyhedral_fan(red_scones; non_redundant=true)

    open(XzCompressorStream, output_file, "w+") do out
      save(out, (graph=G,stratum=red_sfan))
    end
  else
    open(XzCompressorStream, output_file, "w+") do out
      save(out, (graph=G,stratum=simplex(0)|>normal_fan))
    end
  end
end
