using JSON
using Oscar
using CodecXz

function parse_triangulation_to_mrdi(dat_file::String, xz_file::String, outfile::String)
  entire = read(dat_file, String)
  m = match(r"(\[\D*\[.*\]\D*\])\D*(\[\D*\[.*\]\D*\])"s, entire)
  if !isnothing(m)
    points = m.captures[1]
  else
    points = entire
  end
  points = matrix(QQ, Meta.eval(Meta.parse(points)))

  try
    global group = m.captures[2]
    global group = Polymake.to_one_based_indexing(convert(Vector{Vector{Int}}, Meta.parse(group)|>Meta.eval))
    global group = permutation_group(nrows(points), [perm(x) for x in group])
  catch e
    global group = permutation_group(nrows(points), PermGroupElem[])
  end

  triangs = IncidenceMatrix[]
  global i = 1
  open(`xzcat $xz_file`, "r") do io
      while !eof(io)
          line = readline(io)
          m = match(r"{{.*}}", line)
          triang = replace(m.match, "{" => "[")
          triang = replace(triang, "}" => "]")
          triang = convert(Vector{Vector{Int}}, JSON.parse(triang))
          triang = Polymake.to_one_based_indexing(triang)
          push!(triangs, IncidenceMatrix(triang))
          global i += 1
      end
  end
  println("Found $(i-1) triangulations in $xz_file")
  tuple = (points=points, group=group, triangulations=triangs)
  open(XzCompressorStream,outfile,"w+") do out
    save(out, tuple)
  end
end

