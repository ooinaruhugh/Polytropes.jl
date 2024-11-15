module Polytropes

import Oscar

include("Graphs.jl")
include("weighted_digraph_polyhedron.jl")

include("Cayley.jl")

include("Transitivity.jl")
import .Transitivity: transitive_closure, transitive_reduction
export transitive_closure
export transitive_reduction
export transitively_closed_acyclic_graphs

export complete_dag
export opposite_graph
export indegree
export outdegree

export weighted_digraph_polyhedron

end # module Polytropes
