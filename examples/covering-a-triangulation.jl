using Oscar, Polytropes
import Polytropes: subdivision_of_fundamental_polytope

function count_incidences(I::IncidenceMatrix)
  n = size(cells,2)

  return [sum(I[:,i]) for i in 1:n]
end

function find_loneliest_points(I::IncidenceMatrix)
  counts = count_incidences(I)
  mindeg = filter(!iszero, counts) |> minimum 

  return findall(==(mindeg), counts)
end

function find_next_set(I::IncidenceMatrix, v::Int, covered::Vector{Int})
  return filter(eachrow(I) |> enumerate |> collect) do (i, row)
    row[v] && any(x -> (first(x)!=1) || (last(x)!=1), zip(row,covered))
  end |> rand
end

function qd_vertex_cover(I::IncidenceMatrix; totally_randomized=false)
  vc = []
  rows = [1:size(I, 1)...]
  n = size(I,2)
  covered = [0 for _ in 1:n]

  while 0 in covered
    next = find_loneliest_points(I[rows,:]) |> rand
    i, row = find_next_set(I, next, covered)

    setdiff!(rows, [i])
    push!(vc,i)
    covered += row
  end

  return I[unique(vc),:]
end

n = 4
G = complete_dag(n)

sop = subdivision_of_fundamental_polytope(G, [0,0,0,3,1,3])

cells = maximal_cells(IncidenceMatrix, sop)
