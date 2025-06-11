module Polytropes

import Oscar

include("graphs.jl")
include("polyhedra.jl")

include("cayley-trick.jl")

include("Transitivity.jl")
import .Transitivity: transitive_closure, transitive_reduction
export transitive_closure
export transitive_reduction
export transitively_closed_acyclic_graphs

export complete_dag

export root_polytope
export fundamental_polytope

export opposite_graph
export indegree
export outdegree

export weighted_digraph_polyhedron

end # module Polytropes
