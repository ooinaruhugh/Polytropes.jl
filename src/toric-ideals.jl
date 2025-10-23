using Oscar
import Oscar: toric_ideal

function variable_labels(G::Graph{Directed})
  E = edges(G)
  return map(e->"e$(src(e))$(dst(e))", E)
end

function edge_ring(G::Graph{Directed})
  E = variable_labels(G)
  R,x = polynomial_ring(QQ,E)

  return R,x
end

function toric_ideal(G::Graph{Directed})
  R,_ = edge_ring(G)
  A = fundamental_polytope(Matrix,G)[2:end,:]

  return toric_ideal(R,A)
end

function initial_form(f::MPolyRingElem, w)
  R = parent(f)
  maxW = -inf
  ctx = MPolyBuildCtx(R)

  for (c,e) in coefficients_and_exponents(f)
    we = dot(w,e)
    if we > maxW
      finish(ctx)
      maxW = we
      push_term!(ctx,c,e)
    elseif we == maxW
      push_term!(ctx,c,w)
    end
  end

  return finish(ctx)
end

